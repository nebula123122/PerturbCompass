from __future__ import annotations
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
import anndata as ad
import h5py
import numpy as np
import pandas as pd
import yaml
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD
from perturbcompass.runtime import ROOT, root_path, digest, verify_migrated_inputs
from perturbcompass import metadata

def choose_probability_selector(losses: np.ndarray, families: list[str], bootstrap_repetitions: int, seed: int, confidence_threshold: float) -> dict[str, Any]:
    """Choose a clear OOF winner or use bootstrap best-family probabilities as weights."""
    if losses.ndim != 2 or losses.shape[1] != len(families) or (not np.isfinite(losses).all()):
        raise ValueError('loss matrix is invalid')
    means = losses.mean(axis=0)
    winner = int(np.argmin(means))
    rng = np.random.default_rng(seed)
    counts = np.zeros(len(families), dtype=np.int64)
    for start in range(0, bootstrap_repetitions, 1000):
        size = min(1000, bootstrap_repetitions - start)
        rows = rng.integers(0, losses.shape[0], size=(size, losses.shape[0]))
        boot = losses[rows].mean(axis=1)
        counts += np.bincount(np.argmin(boot, axis=1), minlength=len(families))
    probabilities = counts.astype(float) / float(bootstrap_repetitions)
    if probabilities[winner] >= confidence_threshold:
        weights = np.zeros(len(families), dtype=float)
        weights[winner] = 1.0
        mode = 'dominant_single_family'
    else:
        weights = probabilities / probabilities.sum()
        mode = 'ambiguous_probability_model_average'
    return {'mode': mode, 'anchor_winner': families[winner], 'leading_probability': float(probabilities[winner]), 'anchor_mean_losses': {family: float(means[i]) for i, family in enumerate(families)}, 'bootstrap_best_probability': {family: float(probabilities[i]) for i, family in enumerate(families)}, 'weights': {family: float(weights[i]) for i, family in enumerate(families)}, 'weight_vector': weights, 'confidence_threshold': float(confidence_threshold)}

def weighted_prediction(predictions: dict[str, np.ndarray], families: list[str], weights: np.ndarray) -> np.ndarray:
    stack = np.stack([predictions[family] for family in families], axis=0)
    return np.tensordot(weights, stack, axes=(0, 0)).astype(np.float32)

def assign_target_programs(targets: list[str], genes: list[str], labels: np.ndarray, source_effects: np.ndarray) -> tuple[np.ndarray, list[str]]:
    """Assign target-gene programs, with a deterministic source-only fallback."""
    gene_to_index = {gene: i for i, gene in enumerate(genes)}
    program_count = int(labels.max()) + 1
    group = (labels[:, None] == np.arange(program_count)[None, :]).astype(float)
    counts = group.sum(axis=0)
    if np.any(counts == 0):
        raise RuntimeError('empty gene program')
    source_energy = source_effects.astype(float) ** 2 @ group / counts[None, :]
    missing: list[str] = []
    assignments = np.empty(len(targets), dtype=int)
    for index, target in enumerate(targets):
        if target in gene_to_index:
            assignments[index] = int(labels[gene_to_index[target]])
        else:
            missing.append(target)
            assignments[index] = int(np.argmax(source_energy[index]))
    return (assignments, missing)

def target_log_fc(path: Path, genes: list[str], targets: list[str], condition: str, layer_name: str) -> tuple[np.ndarray, dict[str, Any]]:
    """Read only one frozen condition and the frozen common genes from a dense H5AD layer."""
    with h5py.File(path, 'r') as handle:
        obs = handle['obs']
        var = handle['var']
        obs_targets = np.asarray(metadata.read_frame_column(obs, 'target_contrast_gene_name'), dtype=str)
        obs_targets = np.char.upper(np.char.strip(obs_targets))
        obs_conditions = np.asarray(metadata.read_frame_column(obs, 'culture_condition'), dtype=str)
        obs_conditions = np.char.strip(obs_conditions)
        measured = metadata.read_frame_column(var, 'gene_name') if 'gene_name' in var else metadata.read_index(var)
        lookup: dict[str, int] = {}
        for index, value in enumerate(measured):
            token = str(value).strip().upper()
            if token and token not in lookup:
                lookup[token] = index
        missing_genes = [gene for gene in genes if gene not in lookup]
        if missing_genes:
            raise RuntimeError(f'target H5AD lost {len(missing_genes)} frozen genes')
        row_map = {str(target): int(index) for index, target in enumerate(obs_targets) if obs_conditions[index] == condition and target in set(targets)}
        missing_targets = [target for target in targets if target not in row_map]
        if missing_targets:
            raise RuntimeError(f'{condition} lacks frozen targets: {missing_targets[:20]}')
        node = handle['layers'][layer_name]
        if not isinstance(node, h5py.Dataset) or node.ndim != 2:
            raise RuntimeError(f'{layer_name} is not a dense HDF5 matrix')
        row_indices = np.asarray([row_map[target] for target in targets], dtype=np.int64)
        order = np.argsort(row_indices, kind='stable')
        sorted_block = np.asarray(node[row_indices[order], :], dtype=np.float32)
        block = np.empty_like(sorted_block)
        block[order] = sorted_block
        columns = np.asarray([lookup[gene] for gene in genes], dtype=np.int64)
        truth = block[:, columns]
    if not np.isfinite(truth).all():
        raise RuntimeError(f'non-finite {layer_name} values in frozen {condition} matrix')
    return (truth, {'condition': condition, 'targets': len(targets), 'genes': len(genes), 'layer': layer_name, 'mean_abs_log_fc': float(np.mean(np.abs(truth))), 'nonzero_fraction': float(np.mean(truth != 0))})
