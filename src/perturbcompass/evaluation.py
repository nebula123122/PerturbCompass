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
import scipy.sparse as sp
import yaml
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD
from perturbcompass.runtime import ROOT, root_path, digest, verify_migrated_inputs

def bootstrap_regret(selector: np.ndarray, family_losses: np.ndarray, repetitions: int, seed: int) -> tuple[np.ndarray, float, list[float]]:
    actual = float(selector.mean() / family_losses.mean(axis=0).min() - 1.0)
    rng = np.random.default_rng(seed)
    values = np.empty(repetitions, dtype=float)
    for start in range(0, repetitions, 1000):
        stop = min(start + 1000, repetitions)
        rows = rng.integers(0, len(selector), size=(stop - start, len(selector)))
        selected_mean = selector[rows].mean(axis=1)
        family_mean = family_losses[rows].mean(axis=1)
        values[start:stop] = selected_mean / family_mean.min(axis=1) - 1.0
    ci = np.quantile(values, [0.025, 0.975]).tolist()
    return (values, actual, [float(ci[0]), float(ci[1])])

def paired_inference(candidate: np.ndarray, reference: np.ndarray, bootstrap_repetitions: int, bootstrap_seed: int, sign_flip_repetitions: int, sign_seed: int) -> dict[str, Any]:
    rng = np.random.default_rng(bootstrap_seed)
    improvement = np.empty(bootstrap_repetitions, dtype=float)
    for start in range(0, bootstrap_repetitions, 1000):
        stop = min(start + 1000, bootstrap_repetitions)
        rows = rng.integers(0, len(candidate), size=(stop - start, len(candidate)))
        improvement[start:stop] = 1.0 - candidate[rows].mean(axis=1) / reference[rows].mean(axis=1)
    difference = reference - candidate
    observed = float(difference.mean())
    rng = np.random.default_rng(sign_seed)
    exceed = 0
    for start in range(0, sign_flip_repetitions, 2000):
        size = min(2000, sign_flip_repetitions - start)
        signs = rng.choice(np.asarray([-1.0, 1.0]), size=(size, len(difference)))
        exceed += int(np.count_nonzero(np.mean(signs * difference[None, :], axis=1) >= observed))
    ci = np.quantile(improvement, [0.025, 0.975])
    return {'relative_improvement': float(1.0 - candidate.mean() / reference.mean()), 'bootstrap_ci95': [float(ci[0]), float(ci[1])], 'sign_flip_one_sided_p': float((exceed + 1) / (sign_flip_repetitions + 1)), 'targets_improved_fraction': float(np.mean(candidate < reference)), 'paired_mean_loss_reduction': observed}

def holm_adjust(p_values: dict[str, float]) -> dict[str, float]:
    ordered = sorted(p_values, key=p_values.get)
    adjusted: dict[str, float] = {}
    running = 0.0
    total = len(ordered)
    for rank, key in enumerate(ordered):
        running = max(running, min(1.0, (total - rank) * p_values[key]))
        adjusted[key] = running
    return adjusted

def scope_indices(source: np.ndarray, count: int | None) -> list[np.ndarray]:
    if count is None or count >= source.shape[1]:
        return [np.arange(source.shape[1], dtype=int) for _ in range(len(source))]
    return [np.argsort(np.abs(row))[-count:] for row in source]

def per_target_mse(prediction: np.ndarray, truth: np.ndarray, indices: list[np.ndarray]) -> np.ndarray:
    return np.asarray([float(np.mean((prediction[i, selected] - truth[i, selected]) ** 2)) for i, selected in enumerate(indices)], dtype=float)
