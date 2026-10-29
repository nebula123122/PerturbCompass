from __future__ import annotations
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
import anndata as ad
import numpy as np
import pandas as pd
import yaml
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from perturbcompass.runtime import ROOT, root_path, digest, verify_migrated_inputs
from perturbcompass import transfer

def stratified_round_robin_folds(values: np.ndarray, fold_count: int, seed: int) -> np.ndarray:
    """Balance source-effect quintiles across seeded target folds."""
    values = np.asarray(values, dtype=float)
    ranks = pd.Series(values).rank(method='first')
    bin_count = min(5, len(values))
    bins = pd.qcut(ranks, q=bin_count, labels=False, duplicates='drop').to_numpy(dtype=int)
    rng = np.random.default_rng(seed)
    folds = np.empty(len(values), dtype=np.int16)
    offset = 0
    for bin_value in sorted(set(bins.tolist())):
        indices = np.flatnonzero(bins == bin_value)
        indices = rng.permutation(indices)
        folds[indices] = (np.arange(len(indices)) + offset) % fold_count
        offset = (offset + len(indices)) % fold_count
    return folds

def membership(labels: np.ndarray, program_count: int) -> tuple[np.ndarray, np.ndarray]:
    matrix = (labels[:, None] == np.arange(program_count)[None, :]).astype(np.float64)
    counts = matrix.sum(axis=0)
    if np.any(counts == 0):
        raise RuntimeError('empty gene program')
    return (matrix, counts)

def program_features_fast(source: np.ndarray, target_indices: np.ndarray, labels: np.ndarray, target_programs: np.ndarray, task_x: np.ndarray, program_count: int) -> tuple[np.ndarray, np.ndarray]:
    selected = source[target_indices].astype(np.float64, copy=False)
    group, counts = membership(labels, program_count)
    energy = selected * selected @ group
    energy = np.maximum(energy, 1e-12)
    total_energy = np.maximum(np.sum(selected * selected, axis=1, keepdims=True), 1e-12)
    rms = np.sqrt(energy / counts[None, :])
    mean_abs = np.abs(selected) @ group / counts[None, :]
    positive = (selected > 0).astype(np.float64) @ group / counts[None, :]
    energy_share = energy / total_energy
    gene_fraction = np.broadcast_to(counts[None, :] / source.shape[1], energy.shape)
    regulatory = target_programs[target_indices, None] == np.arange(program_count)[None, :]
    repeated_task = np.repeat(task_x[target_indices, None, :], program_count, axis=1)
    one_hot = np.broadcast_to(np.eye(program_count)[None, :, :], (len(target_indices), program_count, program_count))
    features = np.concatenate([rms[:, :, None], mean_abs[:, :, None], positive[:, :, None], energy_share[:, :, None], gene_fraction[:, :, None], regulatory[:, :, None].astype(float), repeated_task, one_hot], axis=2)
    return (features.reshape(-1, features.shape[2]), energy.reshape(-1))

def fit_program_coefficients(source: np.ndarray, target_values: np.ndarray, train: np.ndarray, test: np.ndarray, labels: np.ndarray, target_programs: np.ndarray, task_x: np.ndarray, program_count: int, ridge_alpha: float, bounds: tuple[float, float], centering_strength: float, base_test: np.ndarray | None=None) -> np.ndarray:
    """Predict direct alphas or residual betas using program-level weighted Ridge."""
    train_indices, test_indices = (np.flatnonzero(train), np.flatnonzero(test))
    group, _ = membership(labels, program_count)
    train_source = source[train_indices].astype(np.float64, copy=False)
    numerator = train_source * target_values[train_indices] @ group
    denominator = np.maximum(train_source * train_source @ group, 1e-12)
    oracle = np.clip(numerator / denominator, bounds[0], bounds[1])
    center = np.sum(oracle * denominator, axis=1, keepdims=True) / np.sum(denominator, axis=1, keepdims=True)
    oracle = oracle - centering_strength * center
    x_train, energy_train = program_features_fast(source, train_indices, labels, target_programs, task_x, program_count)
    x_test, energy_test = program_features_fast(source, test_indices, labels, target_programs, task_x, program_count)
    scaler = StandardScaler().fit(x_train)
    model = Ridge(alpha=ridge_alpha).fit(scaler.transform(x_train), oracle.reshape(-1), sample_weight=energy_train / max(float(energy_train.mean()), 1e-12))
    raw = np.asarray(model.predict(scaler.transform(x_test)), dtype=float).reshape(len(test_indices), program_count)
    energy_matrix = energy_test.reshape(len(test_indices), program_count)
    raw_center = np.sum(raw * energy_matrix, axis=1, keepdims=True) / np.sum(energy_matrix, axis=1, keepdims=True)
    raw = raw - centering_strength * raw_center
    if base_test is not None:
        raw = raw + np.asarray(base_test)[:, None]
    return np.clip(raw, bounds[0], bounds[1])

def inner_oof_low_rank(source: np.ndarray, truth: np.ndarray, outer_train: np.ndarray, source_rms: np.ndarray, config: dict[str, Any], outer_fold: int) -> np.ndarray:
    train_indices = np.flatnonzero(outer_train)
    local_folds = stratified_round_robin_folds(source_rms[train_indices], int(config['inner_target_folds']), int(config['inner_fold_seed']) + outer_fold)
    predictions = np.zeros_like(truth, dtype=np.float32)
    for inner_fold in range(int(config['inner_target_folds'])):
        inner_test_indices = train_indices[local_folds == inner_fold]
        inner_test = np.zeros(len(source), dtype=bool)
        inner_test[inner_test_indices] = True
        inner_train = outer_train & ~inner_test
        predictions[inner_test] = transfer.low_rank_prediction(source, truth, inner_train, inner_test, int(config['low_rank_components']), float(config['low_rank_ridge_alpha']), int(config['inner_fold_seed']) + 10 * outer_fold + inner_fold)
    return predictions

def infer_pair(candidate: np.ndarray, reference: np.ndarray, comparator: str, repetitions_bootstrap: int, bootstrap_seed: int, repetitions_sign_flip: int, sign_flip_seed: int) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    rng = np.random.default_rng(bootstrap_seed)
    improvements = np.empty(repetitions_bootstrap, dtype=float)
    bootstrap_rows = []
    for index in range(repetitions_bootstrap):
        selected = rng.integers(0, len(candidate), size=len(candidate))
        improvements[index] = 1.0 - candidate[selected].mean() / reference[selected].mean()
        bootstrap_rows.append({'comparator': comparator, 'bootstrap': index, 'improvement': improvements[index]})
    difference = reference - candidate
    observed = float(difference.mean())
    rng = np.random.default_rng(sign_flip_seed)
    null = np.empty(repetitions_sign_flip, dtype=float)
    for start in range(0, repetitions_sign_flip, 2000):
        stop = min(start + 2000, repetitions_sign_flip)
        signs = rng.choice([-1.0, 1.0], size=(stop - start, len(difference)))
        null[start:stop] = np.mean(signs * difference[None, :], axis=1)
    sign_rows = [{'comparator': comparator, 'permutation': index, 'null_mean_loss_difference': value} for index, value in enumerate(null)]
    ci = np.quantile(improvements, [0.025, 0.975])
    p_value = float((1 + np.count_nonzero(null >= observed)) / (1 + repetitions_sign_flip))
    sd = float(np.std(difference, ddof=1))
    return ({'comparator': comparator, 'relative_improvement': float(1.0 - candidate.mean() / reference.mean()), 'bootstrap_ci95': [float(ci[0]), float(ci[1])], 'sign_flip_one_sided_p': p_value, 'targets_improved_fraction': float(np.mean(candidate < reference)), 'paired_mean_loss_reduction': observed, 'paired_loss_reduction_sd': sd, 'cohen_dz': observed / sd if sd > 0 else None}, bootstrap_rows, sign_rows)
