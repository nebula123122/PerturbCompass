from __future__ import annotations
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import anndata as ad
import h5py
import numpy as np
import pandas as pd
import yaml
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import normalize
from perturbcompass.runtime import ROOT, root_path, digest, verify_migrated_inputs

def relabel_clusters(labels: np.ndarray, centers: np.ndarray) -> np.ndarray:
    ordering = sorted(range(len(centers)), key=lambda index: (float(centers[index, 0]), index))
    mapping = {old: new for new, old in enumerate(ordering)}
    return np.asarray([mapping[int(label)] for label in labels], dtype=int)

def cluster_counts(partition: pd.DataFrame, labels: np.ndarray, cluster_count: int) -> list[dict[str, int]]:
    rows = []
    for cluster in range(cluster_count):
        in_cluster = labels == cluster
        training_anchors = int((partition['role'].eq('anchor').to_numpy() & ~in_cluster).sum())
        test_targets = int((partition['role'].eq('test').to_numpy() & in_cluster).sum())
        rows.append({'cluster': cluster, 'training_anchors_outside_cluster': training_anchors, 'heldout_test_targets': test_targets})
    return rows

def farthest_first_order(features: np.ndarray, names: list[str]) -> list[int]:
    """Deterministic centroid start followed by maximin distance and lexical ties."""
    values = np.asarray(features, dtype=float)
    values = values / np.maximum(np.linalg.norm(values, axis=1, keepdims=True), 1e-12)
    centroid = values.mean(axis=0)
    first_distance = np.sum((values - centroid) ** 2, axis=1)
    first_candidates = np.flatnonzero(np.isclose(first_distance, first_distance.min(), rtol=1e-12, atol=1e-14))
    first = min(first_candidates, key=lambda index: names[int(index)])
    order = [int(first)]
    selected = np.zeros(len(values), dtype=bool)
    selected[first] = True
    minimum = np.sum((values - values[first]) ** 2, axis=1)
    while len(order) < len(values):
        maximum = minimum[~selected].max()
        candidates = np.flatnonzero(~selected & np.isclose(minimum, maximum, rtol=1e-12, atol=1e-14))
        chosen = min(candidates, key=lambda index: names[int(index)])
        order.append(int(chosen))
        selected[chosen] = True
        minimum = np.minimum(minimum, np.sum((values - values[chosen]) ** 2, axis=1))
    return order
