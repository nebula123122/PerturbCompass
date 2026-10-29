from __future__ import annotations
import hashlib
import json
from pathlib import Path
from typing import Any
import anndata as ad
import h5py
import numpy as np
import pandas as pd
import yaml
from perturbcompass.runtime import ROOT, root_path, digest, verify_migrated_inputs

def decode(values: Any) -> np.ndarray:
    array = np.asarray(values)
    if array.dtype.kind in {'S', 'O'}:
        return np.asarray([value.decode('utf-8') if isinstance(value, (bytes, np.bytes_)) else str(value) for value in array], dtype=object)
    return array

def read_frame_column(group: h5py.Group, key: str) -> np.ndarray:
    node = group[key]
    if isinstance(node, h5py.Dataset):
        return decode(node[...])
    if isinstance(node, h5py.Group) and {'codes', 'categories'}.issubset(node.keys()):
        codes = np.asarray(node['codes'][...], dtype=int)
        categories = decode(node['categories'][...])
        result = np.empty(len(codes), dtype=object)
        result[codes < 0] = None
        valid = codes >= 0
        result[valid] = categories[codes[valid]]
        return result
    raise RuntimeError(f'unsupported H5AD dataframe encoding at {group.name}/{key}')

def read_index(group: h5py.Group) -> np.ndarray:
    key = group.attrs.get('_index', '_index')
    if isinstance(key, bytes):
        key = key.decode('utf-8')
    return read_frame_column(group, str(key))

def allocate_anchors(frame: pd.DataFrame, total: int, seed: int) -> set[str]:
    grouped = frame.groupby('stratum', sort=True)
    quotas: dict[str, int] = {}
    remainders: list[tuple[float, str]] = []
    for stratum, group in grouped:
        exact = len(group) * total / len(frame)
        quotas[str(stratum)] = int(np.floor(exact))
        remainders.append((exact - np.floor(exact), str(stratum)))
    remaining = total - sum(quotas.values())
    for _, stratum in sorted(remainders, key=lambda item: (-item[0], item[1]))[:remaining]:
        quotas[stratum] += 1
    rng = np.random.default_rng(seed)
    selected: set[str] = set()
    for stratum, group in grouped:
        values = np.asarray(sorted(group['target'].astype(str)), dtype=object)
        selected.update(map(str, rng.permutation(values)[:quotas[str(stratum)]]))
    if len(selected) != total:
        raise RuntimeError(f'allocated {len(selected)} anchors rather than {total}')
    return selected
