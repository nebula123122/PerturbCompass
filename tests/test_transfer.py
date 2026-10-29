from __future__ import annotations
from perturbcompass import transfer as transfer
from perturbcompass import transfer, programs, evaluation, apma, design as prepare
from pathlib import Path
import numpy as np
import pandas as pd

def test_adamson_label_cleaning_matches_frozen_rules() -> None:
    raw = pd.Series(['SEC61A1_pDS031', '63(mod)_pBA580', '*', 'nan'])
    config = {'truncate_at': '_', 'control_regex': '\\(', 'drop_labels': ['*', 'nan', 'NaN', 'None', '']}
    assert transfer.clean_target_labels(raw, config).tolist() == ['SEC61A1', 'CONTROL', 'DROP', 'DROP']

def test_holm_adjust_is_monotone_in_sorted_p_values() -> None:
    raw = [0.04, 0.001, 0.02]
    adjusted = transfer.holm_adjust(raw)
    ordered = np.argsort(raw)
    assert all((adjusted[ordered[index]] <= adjusted[ordered[index + 1]] for index in range(len(raw) - 1)))
    assert all((0.0 <= value <= 1.0 for value in adjusted))

def test_gene_wise_ridge_returns_test_by_gene_matrix() -> None:
    source = np.arange(30, dtype=float).reshape(6, 5)
    truth = 0.5 * source + 1.0
    train = np.array([True, True, True, True, False, False])
    prediction = transfer.gene_wise_prediction(source, truth, train, ~train, ridge_alpha=0.0)
    assert prediction.shape == (2, 5)
    assert np.allclose(prediction, truth[~train])

def test_program_source_features_preserve_owner_blocks() -> None:
    source = np.arange(24, dtype=float).reshape(3, 8)
    labels = np.repeat(np.arange(2), 4)
    target_programs = np.array([0, 1, 0])
    task_x = np.ones((3, 5))
    x, energy, owner = transfer.program_source_features(source, np.array([0, 2]), labels, target_programs, task_x, 2)
    assert x.shape[0] == 4
    assert energy.shape == (4,)
    assert owner.tolist() == [0, 0, 2, 2]
