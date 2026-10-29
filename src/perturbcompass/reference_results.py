"""Independent, dependency-free checks of retained paper result tables."""
from __future__ import annotations

import csv
import json
import math
from statistics import fmean

from perturbcompass.runtime import ROOT, verify_small_artifacts


def _table(path):
    with (ROOT / path).open(newline='', encoding='utf-8') as stream:
        return list(csv.DictReader(stream))


def audit_reference_results():
    verify_small_artifacts()
    confirmation = _table('results/reference/cd4t_confirmation/heldout_metrics.csv')
    budget = _table('results/reference/anchor_budget/heldout_budget_metrics.csv')
    summary = json.loads((ROOT / 'results/reference/cd4t_confirmation/cd4t_prospective_summary.json').read_text())
    families = ['direct_program', 'scalar', 'residual_program', 'low_rank', 'gene_ridge']
    assert len(confirmation) == 1455 and len(budget) == 5820
    lookup = {(r['condition'], r['target']): r for r in budget if int(r['budget']) == 256}
    max_difference = max(abs(float(row[column]) - float(lookup[row['condition'], row['target']][column]))
                         for row in confirmation for column in [name + '_mse_top200' for name in families + ['selector']])
    assert max_difference <= 1e-9
    results = {}
    for condition in ['Rest', 'Stim8hr', 'Stim48hr']:
        rows = [r for r in confirmation if r['condition'] == condition]
        assert len(rows) == 485 and len({r['target'] for r in rows}) == 485
        means = {family: fmean(float(r[family + '_mse_top200']) for r in rows) for family in families}
        loss = fmean(float(r['selector_mse_top200']) for r in rows)
        regret = loss / min(means.values()) - 1
        recorded = summary['conditions'][condition]['heldout']['selector_regret']
        assert math.isclose(regret, recorded, rel_tol=1e-10, abs_tol=1e-12)
        small = [r for r in budget if r['condition'] == condition and int(r['budget']) == 32]
        increase = fmean(float(r['selector_mse_top200']) for r in small) / loss - 1
        results[condition] = {'heldout_targets': 485, 'selector_regret': regret, 'anchor_32_vs_256_relative_mse_increase': increase}
    return {'status': 'REFERENCE_RESULTS_RECONCILED', 'conditions': results, 'full_budget_max_absolute_difference': max_difference}
