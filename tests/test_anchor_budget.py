from perturbcompass.analyses import anchor_budget as budget
from perturbcompass import transfer, programs, evaluation, apma, design as prepare
from pathlib import Path
import numpy as np

def test_farthest_first_is_deterministic_permutation_and_nested():
    features = np.asarray([[0.0, 0.0], [1.0, 0.0], [-1.0, 0.0], [0.0, 2.0]])
    names = ['D', 'B', 'A', 'C']
    first = prepare.farthest_first_order(features, names)
    second = prepare.farthest_first_order(features, names)
    assert first == second
    assert sorted(first) == list(range(4))
    assert set(first[:2]).issubset(set(first[:3]))

def test_relative_increase_bootstrap_contains_point_estimate():
    reference = np.linspace(1.0, 2.0, 50)
    candidate = reference * 1.02
    actual, interval = budget.relative_increase_bootstrap(candidate, reference, 4000, 1)
    assert np.isclose(actual, 0.02)
    assert interval[0] <= actual <= interval[1]

def test_relative_increase_zero_for_identical_losses():
    values = np.linspace(0.5, 1.5, 30)
    actual, interval = budget.relative_increase_bootstrap(values, values, 2000, 2)
    assert np.isclose(actual, 0.0)
    assert np.allclose(interval, [0.0, 0.0])

def test_canonical_anchor_indices_ignore_selection_rank():
    lookup = {'A': 0, 'B': 1, 'C': 2, 'D': 3}
    observed = budget.canonical_anchor_indices(['D', 'B', 'C'], lookup)
    assert observed.tolist() == [1, 2, 3]
