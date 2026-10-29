from perturbcompass.analyses import source_ood as ood
from perturbcompass import transfer, programs, evaluation, apma, design as prepare
from pathlib import Path
import numpy as np

def test_stratified_bootstrap_preserves_cluster_sizes_and_actual_regret():
    selector = np.asarray([1.0, 1.0, 2.0, 2.0, 3.0, 3.0])
    families = np.column_stack([selector + 0.2, selector + 0.4])
    clusters = np.asarray([0, 0, 1, 1, 2, 2])
    values, actual, interval = ood.stratified_bootstrap_regret(selector, families, clusters, 2000, 1)
    expected = selector.mean() / (selector + 0.2).mean() - 1.0
    assert np.isclose(actual, expected)
    assert len(values) == 2000
    assert interval[0] <= actual <= interval[1]

def test_cluster_relabel_is_deterministic():
    from perturbcompass import design as prepare
    labels = np.asarray([2, 0, 1, 2])
    centers = np.asarray([[2.0, 0.0], [-1.0, 0.0], [0.5, 0.0]])
    observed = prepare.relabel_clusters(labels, centers)
    assert observed.tolist() == [1, 2, 0, 1]

def test_selector_probability_weights_remain_convex():
    losses = np.asarray([[1.0, 1.1], [1.2, 1.0], [1.0, 1.2], [1.2, 1.0]])
    result = apma.choose_probability_selector(losses, ['a', 'b'], 2000, 2, 0.8)
    assert np.isclose(sum(result['weights'].values()), 1.0)
    assert all((value >= 0 for value in result['weights'].values()))

def test_unmeasured_target_program_falls_back_without_outcome_access():
    targets = ['A', 'MISSING']
    genes = ['A', 'B', 'C', 'D']
    gene_programs = np.asarray([0, 0, 1, 1])
    source = np.asarray([[2.0, 1.0, 0.0, 0.0], [0.0, 0.0, 2.0, 1.0]])
    assigned, fallback = apma.assign_target_programs(targets, genes, gene_programs, source)
    assert assigned.tolist() == [0, 1]
    assert fallback == ['MISSING']

def test_cluster_counts_use_disjoint_roles():
    import pandas as pd
    from perturbcompass import design as prepare
    partition = pd.DataFrame({'role': ['anchor', 'anchor', 'test', 'test', 'test']})
    rows = prepare.cluster_counts(partition, np.asarray([0, 1, 0, 1, 1]), 2)
    assert rows[0]['training_anchors_outside_cluster'] == 1
    assert rows[0]['heldout_test_targets'] == 1
    assert rows[1]['training_anchors_outside_cluster'] == 1
    assert rows[1]['heldout_test_targets'] == 2
