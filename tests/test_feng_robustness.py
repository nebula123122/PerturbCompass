from perturbcompass.analyses import feng_robustness as robustness
from perturbcompass import transfer, programs, evaluation, apma, design as prepare
from pathlib import Path
import numpy as np

def test_relative_improvement_direction():
    candidate = np.asarray([1.0, 2.0])
    comparator = np.asarray([2.0, 4.0])
    assert np.isclose(robustness.relative_improvement(candidate, comparator), 0.5)

def test_largest_influence_removal_is_deterministic():
    candidate = np.asarray([1.0, 1.0, 9.0, 1.0])
    comparator = np.asarray([2.0, 2.0, 1.0, 2.0])
    rows = robustness.leave_largest_out(candidate, comparator, [1])
    assert rows[0]['remaining'] == 3
    assert rows[0]['relative_improvement'] > 0

def test_trimmed_improvement_ignores_symmetric_extremes():
    candidate = np.asarray([0.0, 1.0, 1.0, 1.0, 10.0])
    comparator = np.asarray([0.0, 2.0, 2.0, 2.0, 20.0])
    value = robustness.trimmed_relative_improvement(candidate, comparator, 0.2)
    assert np.isclose(value, 0.5)
