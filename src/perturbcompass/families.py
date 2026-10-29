from __future__ import annotations
import argparse
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
from perturbcompass.runtime import ROOT, root_path, digest, verify_migrated_inputs
from perturbcompass import transfer, programs

def family_predictions(source: np.ndarray, truth: np.ndarray, local_train: np.ndarray, local_test: np.ndarray, source_rms: np.ndarray, labels: np.ndarray, target_programs: np.ndarray, task_x: np.ndarray, signatures: list[np.ndarray], model: dict[str, Any], fold: int) -> dict[str, np.ndarray]:
    low_rank = transfer.low_rank_prediction(source, truth, local_train, local_test, int(model['low_rank_components']), float(model['low_rank_ridge_alpha']), int(model['inner_fold_seed']) + fold)
    gene_ridge = transfer.gene_wise_prediction(source, truth, local_train, local_test, float(model['low_rank_ridge_alpha']))
    lower, upper = map(float, model['direct_program_alpha_bounds'])
    oracles, curvature = transfer.scalar_oracles(source, truth, signatures, lower, upper)
    scalar = np.clip(transfer.fit_weighted_ridge(task_x, oracles, curvature, local_train, local_test, float(model['scalar_ridge_alpha'])), lower, upper)
    direct_config = {'program_count': int(model['program_count']), 'alpha_bounds': model['direct_program_alpha_bounds'], 'program_ridge_alpha': float(model['program_ridge_alpha']), 'centering_strength': float(model['direct_program_centering_strength']), 'delta_scale': 1.0}
    direct_alpha = transfer.fit_program_alphas(source, truth, local_train, local_test, labels, target_programs, task_x, scalar, direct_config)
    inner_prediction = programs.inner_oof_low_rank(source, truth, local_train, source_rms, model, outer_fold=100 + fold)
    residual_truth = truth - inner_prediction
    residual_beta = programs.fit_program_coefficients(source, residual_truth, local_train, local_test, labels, target_programs, task_x, int(model['program_count']), float(model['program_ridge_alpha']), tuple(map(float, model['residual_program_beta_bounds'])), float(model['residual_centering_strength']))
    test_indices = np.flatnonzero(local_test)
    return {'low_rank': low_rank, 'gene_ridge': gene_ridge, 'scalar': scalar[:, None] * source[test_indices], 'direct_program': direct_alpha[:, labels] * source[test_indices], 'residual_program': low_rank + residual_beta[:, labels] * source[test_indices]}
