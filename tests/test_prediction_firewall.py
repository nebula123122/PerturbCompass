"""Held-out truth must not affect fitted cross-context predictions."""
import numpy as np

from perturbcompass import evaluation, families, transfer


def test_all_five_families_ignore_heldout_truth():
    rng = np.random.default_rng(981)
    source = rng.normal(size=(30, 48)).astype(np.float32)
    truth = (0.3 * source + rng.normal(0, 0.02, source.shape)).astype(np.float32)
    train = np.arange(30) < 22
    labels = np.repeat(np.arange(4), 12)
    owner = np.arange(30) % 4
    rms = np.sqrt(np.mean(source * source, axis=1))
    features = transfer.task_features(source, np.full(30, 80), np.ones(30))
    scopes = evaluation.scope_indices(source, 20)
    model = {
        'program_count': 4, 'inner_target_folds': 4, 'inner_fold_seed': 303002,
        'low_rank_components': 4, 'low_rank_ridge_alpha': 10., 'scalar_ridge_alpha': 100.,
        'program_ridge_alpha': 10., 'direct_program_alpha_bounds': [0., 1.5],
        'direct_program_centering_strength': .5, 'residual_program_beta_bounds': [-.75, .75],
        'residual_centering_strength': 1.,
    }
    first = families.family_predictions(source, truth, train, ~train, rms, labels, owner, features, scopes, model, 2)
    truth[~train] = 1e5 * rng.normal(size=truth[~train].shape)
    second = families.family_predictions(source, truth, train, ~train, rms, labels, owner, features, scopes, model, 2)
    for name in first:
        np.testing.assert_array_equal(first[name], second[name])
