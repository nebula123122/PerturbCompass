from __future__ import annotations
import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from perturbcompass.runtime import ROOT, root_path, digest, verify_migrated_inputs

def parse_source_target(value: object) -> str:
    parts = str(value).split('_')
    if len(parts) < 2 or not parts[0].isdigit() or (not parts[1]):
        raise ValueError(f'unrecognized source perturbation index: {value}')
    return parts[1].upper()

def gene_lookup(data: ad.AnnData) -> dict[str, int]:
    values = data.var['gene_name'].astype(str) if 'gene_name' in data.var else pd.Index(data.var_names).astype(str)
    result: dict[str, int] = {}
    for index, name in enumerate(values):
        upper = str(name).upper()
        if upper and upper != 'NAN' and (upper not in result):
            result[upper] = index
    return result

def clean_target_labels(raw: pd.Series, config: dict[str, Any]) -> np.ndarray:
    values = raw.astype(str)
    labels = values.str.split(str(config['truncate_at']), n=1).str[0]
    control = values.str.contains(str(config['control_regex']), regex=True)
    labels = labels.mask(control, 'CONTROL')
    labels = labels.mask(values.isin([str(value) for value in config['drop_labels']]), 'DROP')
    return labels.str.upper().to_numpy(dtype=object)

def dense_slice(matrix: Any, columns: np.ndarray) -> np.ndarray:
    selected = matrix[:, columns]
    if sp.issparse(selected):
        selected = selected.toarray()
    return np.asarray(selected, dtype=np.float32)

def aggregate_source(data: ad.AnnData, columns: np.ndarray, minimum_cells: int) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    labels = np.asarray([parse_source_target(value) for value in data.obs_names], dtype=object)
    cells = data.obs['num_cells_filtered'].to_numpy(dtype=float)
    matrix = dense_slice(data.X, columns)
    keep = (cells >= minimum_cells) & np.isfinite(cells) & np.isfinite(matrix).all(axis=1)
    labels, cells, matrix = (labels[keep], cells[keep], matrix[keep])
    targets = sorted(set(labels))
    means = np.empty((len(targets), matrix.shape[1]), dtype=np.float32)
    total_cells = np.empty(len(targets), dtype=float)
    constructs = np.empty(len(targets), dtype=np.int64)
    for index, target in enumerate(targets):
        selected = labels == target
        weights = cells[selected]
        means[index] = np.average(matrix[selected], axis=0, weights=weights).astype(np.float32)
        total_cells[index] = float(weights.sum())
        constructs[index] = int(selected.sum())
    return (targets, means, total_cells, constructs)

def target_statistics(target_data: ad.AnnData, genes: list[str], target_config: dict[str, Any]) -> dict[str, Any]:
    lookup = gene_lookup(target_data)
    columns = np.asarray([lookup[gene] for gene in genes], dtype=np.int64)
    labels = clean_target_labels(target_data.obs[str(target_config['label_column'])], target_config)
    ngenes = target_data.obs['ngenes'].to_numpy(dtype=float)
    ncounts = target_data.obs['ncounts'].to_numpy(dtype=float)
    keep = (labels != 'DROP') & (ngenes >= int(target_config['minimum_genes_per_cell'])) & np.isfinite(ncounts) & (ncounts > 0)
    matrix = target_data.X[keep][:, columns]
    if not sp.issparse(matrix):
        matrix = sp.csr_matrix(matrix)
    matrix = matrix.tocsr().astype(np.float32)
    scale = float(target_config['normalize_total_target_sum']) / ncounts[keep]
    matrix = matrix.multiply(scale[:, None]).tocsr()
    matrix.data = np.log1p(matrix.data).astype(np.float32)
    kept_labels = labels[keep]
    control_rows = np.flatnonzero(kept_labels == 'CONTROL')
    if len(control_rows) < 100:
        raise RuntimeError('too few target control cells after QC')
    fold_count = int(target_config['control_folds'])
    rng = np.random.default_rng(int(target_config['control_fold_seed']))
    shuffled = rng.permutation(control_rows)
    control_fold = np.empty(len(control_rows), dtype=np.int16)
    control_fold[np.arange(len(control_rows))] = np.arange(len(control_rows)) % fold_count
    row_to_fold = dict(zip(shuffled.tolist(), control_fold.tolist(), strict=True))
    control_test_means = []
    control_train_means = []
    control_fold_sizes = []
    for fold in range(fold_count):
        test_rows = np.asarray([row for row in control_rows if row_to_fold[int(row)] == fold], dtype=np.int64)
        train_rows = np.asarray([row for row in control_rows if row_to_fold[int(row)] != fold], dtype=np.int64)
        control_test_means.append(np.asarray(matrix[test_rows].mean(axis=0)).ravel().astype(np.float32))
        control_train_means.append(np.asarray(matrix[train_rows].mean(axis=0)).ravel().astype(np.float32))
        control_fold_sizes.append({'fold': fold, 'train': len(train_rows), 'test': len(test_rows)})
    target_means: dict[str, np.ndarray] = {}
    target_cells: dict[str, int] = {}
    minimum = int(target_config['minimum_target_cells_after_qc'])
    for target in sorted(set(kept_labels) - {'CONTROL', 'DROP'}):
        rows = np.flatnonzero(kept_labels == target)
        if len(rows) >= minimum:
            target_means[target] = np.asarray(matrix[rows].mean(axis=0)).ravel().astype(np.float32)
            target_cells[target] = int(len(rows))
    return {'genes': genes, 'target_means': target_means, 'target_cells': target_cells, 'control_train_means': control_train_means, 'control_test_means': control_test_means, 'control_fold_sizes': control_fold_sizes, 'qc_cells': int(keep.sum()), 'control_cells': int(len(control_rows))}

def task_features(source: np.ndarray, cells: np.ndarray, constructs: np.ndarray) -> np.ndarray:
    rms = np.sqrt(np.maximum(np.mean(source * source, axis=1), 1e-12))
    return np.column_stack([rms, np.mean(np.abs(source), axis=1), np.mean(source > 0, axis=1), np.log1p(cells), np.log1p(constructs)])

def fit_weighted_ridge(x: np.ndarray, y: np.ndarray, weight: np.ndarray, train: np.ndarray, test: np.ndarray, alpha: float) -> np.ndarray:
    scaler = StandardScaler().fit(x[train])
    normalized_weight = weight[train] / max(float(weight[train].mean()), 1e-12)
    model = Ridge(alpha=alpha).fit(scaler.transform(x[train]), y[train], sample_weight=normalized_weight)
    return np.asarray(model.predict(scaler.transform(x[test])), dtype=float)

def scalar_oracles(source: np.ndarray, truth: np.ndarray, signatures: list[np.ndarray], lower: float, upper: float) -> tuple[np.ndarray, np.ndarray]:
    values = np.empty(len(source), dtype=float)
    curvature = np.empty(len(source), dtype=float)
    for index, signature in enumerate(signatures):
        s, t = (source[index, signature], truth[index, signature])
        denominator = max(float(np.dot(s, s)), 1e-12)
        values[index] = np.clip(np.dot(s, t) / denominator, lower, upper)
        curvature[index] = max(float(np.mean(s * s)), 1e-12)
    return (values, curvature)

def program_source_features(source: np.ndarray, target_indices: np.ndarray, labels: np.ndarray, target_programs: np.ndarray, task_x: np.ndarray, program_count: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows, energies, owners = ([], [], [])
    for target_index in target_indices:
        total_energy = max(float(np.dot(source[target_index], source[target_index])), 1e-12)
        for program in range(program_count):
            mask = labels == program
            s = source[target_index, mask]
            energy = max(float(np.dot(s, s)), 1e-12)
            row = [np.sqrt(max(float(np.mean(s * s)), 1e-12)), float(np.mean(np.abs(s))), float(np.mean(s > 0)), energy / total_energy, float(np.mean(mask)), float(program == target_programs[target_index]), *task_x[target_index].tolist(), *[float(number == program) for number in range(program_count)]]
            rows.append(row)
            energies.append(energy)
            owners.append(target_index)
    return (np.asarray(rows, dtype=float), np.asarray(energies, dtype=float), np.asarray(owners, dtype=int))

def fit_program_alphas(source: np.ndarray, truth_train: np.ndarray, train: np.ndarray, test: np.ndarray, labels: np.ndarray, target_programs: np.ndarray, task_x: np.ndarray, scalar_prediction: np.ndarray, config: dict[str, Any]) -> np.ndarray:
    program_count = int(config['program_count'])
    lower, upper = map(float, config['alpha_bounds'])
    train_indices, test_indices = (np.flatnonzero(train), np.flatnonzero(test))
    x_train, energy_train, owner_train = program_source_features(source, train_indices, labels, target_programs, task_x, program_count)
    raw_oracle = np.empty(len(x_train), dtype=float)
    position = 0
    for target_index in train_indices:
        local = []
        for program in range(program_count):
            mask = labels == program
            s, t = (source[target_index, mask], truth_train[target_index, mask])
            value = np.clip(np.dot(s, t) / max(float(np.dot(s, s)), 1e-12), lower, upper)
            local.append(float(value))
        local_array = np.asarray(local)
        local_energy = energy_train[position:position + program_count]
        weighted_mean = float(np.average(local_array, weights=local_energy))
        raw_oracle[position:position + program_count] = local_array - weighted_mean
        position += program_count
    x_test, energy_test, owner_test = program_source_features(source, test_indices, labels, target_programs, task_x, program_count)
    scaler = StandardScaler().fit(x_train)
    model = Ridge(alpha=float(config['program_ridge_alpha'])).fit(scaler.transform(x_train), raw_oracle, sample_weight=energy_train / max(float(energy_train.mean()), 1e-12))
    raw = np.asarray(model.predict(scaler.transform(x_test)), dtype=float)
    result = np.empty((len(test_indices), program_count), dtype=float)
    for local_index, target_index in enumerate(test_indices):
        selected = owner_test == target_index
        delta = raw[selected]
        weights = energy_test[selected]
        delta = delta - float(config['centering_strength']) * float(np.average(delta, weights=weights))
        result[local_index] = np.clip(scalar_prediction[local_index] + float(config['delta_scale']) * delta, lower, upper)
    return result

def global_alpha(source: np.ndarray, truth: np.ndarray, signatures: list[np.ndarray], train: np.ndarray, lower: float, upper: float) -> float:
    numerator, denominator = (0.0, 0.0)
    for index in np.flatnonzero(train):
        s, t = (source[index, signatures[index]], truth[index, signatures[index]])
        numerator += float(np.mean(s * t))
        denominator += float(np.mean(s * s))
    return float(np.clip(numerator / max(denominator, 1e-12), lower, upper))

def gene_wise_prediction(source: np.ndarray, truth: np.ndarray, train: np.ndarray, test: np.ndarray, ridge_alpha: float) -> np.ndarray:
    s, t = (source[train].astype(float), truth[train].astype(float))
    s_mean, t_mean = (s.mean(axis=0), t.mean(axis=0))
    centered_s, centered_t = (s - s_mean, t - t_mean)
    slope = np.sum(centered_s * centered_t, axis=0) / (np.sum(centered_s * centered_s, axis=0) + ridge_alpha)
    intercept = t_mean - slope * s_mean
    return source[test] * slope + intercept

def low_rank_prediction(source: np.ndarray, truth: np.ndarray, train: np.ndarray, test: np.ndarray, components: int, ridge_alpha: float, seed: int) -> np.ndarray:
    rank = max(2, min(components, int(train.sum()) - 1, source.shape[1] - 1))
    source_svd = TruncatedSVD(n_components=rank, random_state=seed).fit(source[train])
    target_svd = TruncatedSVD(n_components=rank, random_state=seed + 1).fit(truth[train])
    x_train = source_svd.transform(source[train])
    x_test = source_svd.transform(source[test])
    y_train = target_svd.transform(truth[train])
    scaler = StandardScaler().fit(x_train)
    model = Ridge(alpha=ridge_alpha).fit(scaler.transform(x_train), y_train)
    return model.predict(scaler.transform(x_test)) @ target_svd.components_

def mse(prediction: np.ndarray, truth: np.ndarray, selected: np.ndarray) -> float:
    return float(np.mean((prediction[selected] - truth[selected]) ** 2))

def holm_adjust(p_values: list[float]) -> list[float]:
    count = len(p_values)
    order = np.argsort(p_values)
    adjusted = np.empty(count, dtype=float)
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (count - rank) * p_values[index]))
        adjusted[index] = running
    return adjusted.tolist()

def infer_pair(candidate: np.ndarray, reference: np.ndarray, source_name: str, comparator: str, config: dict[str, Any], seed_offset: int) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    evaluation = config['evaluation']
    bootstrap_rows = []
    rng = np.random.default_rng(int(evaluation['bootstrap_seed']) + seed_offset)
    improvements = np.empty(int(evaluation['bootstrap_repetitions']), dtype=float)
    for index in range(len(improvements)):
        selected = rng.integers(0, len(candidate), size=len(candidate))
        improvements[index] = 1.0 - candidate[selected].mean() / reference[selected].mean()
        bootstrap_rows.append({'source': source_name, 'comparator': comparator, 'bootstrap': index, 'improvement': improvements[index]})
    difference = reference - candidate
    observed = float(difference.mean())
    repetitions = int(evaluation['sign_flip_repetitions'])
    rng = np.random.default_rng(int(evaluation['sign_flip_seed']) + seed_offset)
    null = np.empty(repetitions, dtype=float)
    for start in range(0, repetitions, 2000):
        stop = min(start + 2000, repetitions)
        signs = rng.choice([-1.0, 1.0], size=(stop - start, len(difference)))
        null[start:stop] = np.mean(signs * difference[None, :], axis=1)
    sign_rows = [{'source': source_name, 'comparator': comparator, 'permutation': index, 'null_mean_loss_difference': value} for index, value in enumerate(null)]
    ci = np.quantile(improvements, [0.025, 0.975])
    p_value = float((1 + np.count_nonzero(null >= observed)) / (1 + repetitions))
    return ({'source': source_name, 'comparator': comparator, 'improvement': float(1.0 - candidate.mean() / reference.mean()), 'ci95': [float(ci[0]), float(ci[1])], 'sign_flip_p': p_value, 'targets_improved_fraction': float(np.mean(candidate < reference))}, bootstrap_rows, sign_rows)
