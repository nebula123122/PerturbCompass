from __future__ import annotations
from perturbcompass import programs as programs
from perturbcompass import transfer, programs, evaluation, apma, design as prepare
from pathlib import Path
import numpy as np
import pandas as pd

def test_stratified_round_robin_folds_are_deterministic_and_balanced() -> None:
    values = np.linspace(0.01, 2.0, 103)
    first = programs.stratified_round_robin_folds(values, 5, 303001)
    second = programs.stratified_round_robin_folds(values, 5, 303001)
    assert np.array_equal(first, second)
    counts = np.bincount(first, minlength=5)
    assert counts.max() - counts.min() <= 1
    ranks = pd.Series(values).rank(method='first')
    quintiles = pd.qcut(ranks, q=5, labels=False, duplicates='drop').to_numpy(dtype=int)
    for quintile in range(5):
        local = np.bincount(first[quintiles == quintile], minlength=5)
        assert local.max() - local.min() <= 1

def test_fast_program_features_match_reference() -> None:
    rng = np.random.default_rng(17)
    source = rng.normal(size=(7, 24))
    labels = np.repeat(np.arange(4), 6)
    targets = np.arange(7) % 4
    task_x = rng.normal(size=(7, 5))
    selected = np.array([0, 2, 5, 6])
    fast_x, fast_energy = programs.program_features_fast(source, selected, labels, targets, task_x, 4)
    reference_x, reference_energy, owners = programs.transfer.program_source_features(source, selected, labels, targets, task_x, 4)
    assert np.array_equal(owners, np.repeat(selected, 4))
    assert np.allclose(fast_x, reference_x, rtol=1e-12, atol=1e-12)
    assert np.allclose(fast_energy, reference_energy, rtol=1e-12, atol=1e-12)

def test_residual_coefficients_are_finite_bounded_and_test_only() -> None:
    rng = np.random.default_rng(23)
    source = rng.normal(size=(12, 32))
    residual = 0.12 * source + rng.normal(scale=0.03, size=source.shape)
    train = np.array([True] * 9 + [False] * 3)
    labels = np.repeat(np.arange(4), 8)
    target_programs = np.arange(12) % 4
    task_x = rng.normal(size=(12, 5))
    beta = programs.fit_program_coefficients(source, residual, train, ~train, labels, target_programs, task_x, 4, 10.0, (-0.75, 0.75), 1.0)
    assert beta.shape == (3, 4)
    assert np.isfinite(beta).all()
    assert np.all(beta >= -0.75) and np.all(beta <= 0.75)

def test_inner_oof_predictions_never_fill_outer_test_rows() -> None:
    rng = np.random.default_rng(29)
    source = rng.normal(size=(15, 40))
    truth = rng.normal(size=(15, 40))
    outer_train = np.array([True] * 11 + [False] * 4)
    config = {'inner_target_folds': 4, 'inner_fold_seed': 303002, 'low_rank_components': 4, 'low_rank_ridge_alpha': 1.0}
    source_rms = np.sqrt(np.mean(source * source, axis=1))
    prediction = programs.inner_oof_low_rank(source, truth, outer_train, source_rms, config, outer_fold=0)
    assert prediction.shape == truth.shape
    assert np.any(prediction[outer_train] != 0)
    assert np.all(prediction[~outer_train] == 0)
