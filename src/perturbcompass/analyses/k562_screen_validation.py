from __future__ import annotations
import argparse
import hashlib
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
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from perturbcompass.runtime import ROOT, root_path, digest, verify_migrated_inputs
from perturbcompass.programs import stratified_round_robin_folds, membership, program_features_fast, fit_program_coefficients, inner_oof_low_rank, infer_pair
from perturbcompass import transfer, programs, families, evaluation, apma, metadata
DEFAULT_CONFIG = ROOT / 'configs/k562_screen_validation.yaml'

def verify_freeze(config_path, config):
    verify_migrated_inputs(__file__, config_path, config)
    return json.loads(root_path(config['outputs']['pre_run_manifest']).read_text(encoding='utf-8'))

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=Path, default=DEFAULT_CONFIG)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text(encoding='utf-8'))
    verify_freeze(config_path, config)
    model = config['model']
    source_data = ad.read_h5ad(root_path(config['inputs']['source_h5ad']))
    target_data = ad.read_h5ad(root_path(config['inputs']['target_h5ad']))
    source_lookup, target_lookup = (transfer.gene_lookup(source_data), transfer.gene_lookup(target_data))
    genes = sorted(set(source_lookup) & set(target_lookup))
    source_columns = np.asarray([source_lookup[gene] for gene in genes], dtype=np.int64)
    target_columns = np.asarray([target_lookup[gene] for gene in genes], dtype=np.int64)
    source_targets, source_means, source_cells, source_constructs = transfer.aggregate_source(source_data, source_columns, int(config['eligibility']['minimum_cells_per_construct']))
    target_targets, target_means, _, _ = transfer.aggregate_source(target_data, target_columns, int(config['eligibility']['minimum_cells_per_construct']))
    source_map = {target: index for index, target in enumerate(source_targets)}
    target_map = {target: index for index, target in enumerate(target_targets)}
    eligible = sorted(set(source_targets) & set(target_targets) & set(genes))
    if len(eligible) < int(config['eligibility']['minimum_eligible_targets']):
        raise RuntimeError(f'only {len(eligible)} eligible targets')
    source = np.stack([source_means[source_map[target]] for target in eligible]).astype(np.float32)
    truth = np.stack([target_means[target_map[target]] for target in eligible]).astype(np.float32)
    cells = np.asarray([source_cells[source_map[target]] for target in eligible], dtype=float)
    constructs = np.asarray([source_constructs[source_map[target]] for target in eligible], dtype=int)
    source_rms = np.sqrt(np.maximum(np.mean(source * source, axis=1), 1e-12))
    folds = stratified_round_robin_folds(source_rms, int(model['outer_target_folds']), int(model['target_fold_seed']))
    embedding = TruncatedSVD(n_components=min(32, len(eligible) - 1, len(genes) - 1), random_state=int(model['program_seed'])).fit_transform(source.T)
    labels = KMeans(n_clusters=int(model['program_count']), random_state=int(model['program_seed']), n_init=20).fit_predict(embedding)
    gene_to_index = {gene: index for index, gene in enumerate(genes)}
    target_programs = np.asarray([labels[gene_to_index[target]] for target in eligible], dtype=int)
    task_x = transfer.task_features(source, cells, constructs)
    signatures = [np.argsort(np.abs(source[index]))[-int(model['top_signature_genes']):] for index in range(len(eligible))]
    random_rng = np.random.default_rng(int(model['random_program_seed']))
    random_labels = [random_rng.permutation(labels) for _ in range(int(model['random_program_partitions']))]
    random_losses = {index: [] for index in range(len(random_labels))}
    rows: list[dict[str, Any]] = []
    lower, upper = map(float, model['direct_program_alpha_bounds'])
    for outer in range(int(model['outer_target_folds'])):
        train, test = (folds != outer, folds == outer)
        low_rank = transfer.low_rank_prediction(source, truth, train, test, int(model['low_rank_components']), float(model['low_rank_ridge_alpha']), int(model['target_fold_seed']) + outer)
        inner_prediction = inner_oof_low_rank(source, truth, train, source_rms, model, outer)
        residual_truth = truth - inner_prediction
        residual_beta = fit_program_coefficients(source, residual_truth, train, test, labels, target_programs, task_x, int(model['program_count']), float(model['program_ridge_alpha']), tuple(map(float, model['residual_program_beta_bounds'])), float(model['residual_centering_strength']))
        oracles, curvature = transfer.scalar_oracles(source, truth, signatures, lower, upper)
        scalar = np.clip(transfer.fit_weighted_ridge(task_x, oracles, curvature, train, test, float(model['scalar_ridge_alpha'])), lower, upper)
        direct_config = {'program_count': int(model['program_count']), 'alpha_bounds': model['direct_program_alpha_bounds'], 'program_ridge_alpha': float(model['program_ridge_alpha']), 'centering_strength': float(model['direct_program_centering_strength']), 'delta_scale': 1.0}
        direct_alpha = transfer.fit_program_alphas(source, truth, train, test, labels, target_programs, task_x, scalar, direct_config)
        gene_ridge = transfer.gene_wise_prediction(source, truth, train, test, float(model['low_rank_ridge_alpha']))
        random_beta = []
        for partition in random_labels:
            random_target_programs = np.asarray([partition[gene_to_index[target]] for target in eligible], dtype=int)
            random_beta.append(fit_program_coefficients(source, residual_truth, train, test, partition, random_target_programs, task_x, int(model['program_count']), float(model['program_ridge_alpha']), tuple(map(float, model['residual_program_beta_bounds'])), float(model['residual_centering_strength'])))
        test_indices = np.flatnonzero(test)
        for local_index, target_index in enumerate(test_indices):
            source_row, target_truth = (source[target_index], truth[target_index])
            predictions = {'low_rank': low_rank[local_index], 'gene_ridge': gene_ridge[local_index], 'direct_program': direct_alpha[local_index][labels] * source_row, 'scalar': scalar[local_index] * source_row, 'residual_program': low_rank[local_index] + residual_beta[local_index][labels] * source_row}
            ranked = np.argsort(np.abs(source_row))
            selections = {'top100': ranked[-min(100, len(ranked)):], 'top200': signatures[target_index], 'top500': ranked[-min(500, len(ranked)):], 'all': np.arange(len(genes))}
            row: dict[str, Any] = {'target': eligible[target_index], 'outer_fold': outer, 'source_cells': cells[target_index], 'source_constructs': constructs[target_index], 'source_effect_rms': source_rms[target_index], 'residual_beta_mean': float(residual_beta[local_index].mean()), 'residual_beta_sd': float(residual_beta[local_index].std())}
            for selection_name, selected in selections.items():
                for method_name, prediction in predictions.items():
                    row[f'{method_name}_mse_{selection_name}'] = transfer.mse(prediction, target_truth, selected)
            rows.append(row)
            selected = selections['top200']
            for partition_index, beta in enumerate(random_beta):
                prediction = low_rank[local_index] + beta[local_index][random_labels[partition_index]] * source_row
                random_losses[partition_index].append(transfer.mse(prediction, target_truth, selected))
        print(json.dumps({'stage': 'v30_outer_fold', 'fold': outer, 'test_targets': int(test.sum())}), flush=True)
    records = pd.DataFrame(rows).sort_values('target').reset_index(drop=True)
    random_records = pd.DataFrame([{'partition': index, 'mean_mse_top200': float(np.mean(losses))} for index, losses in random_losses.items()])
    candidate = records['residual_program_mse_top200'].to_numpy(float)
    comparisons, bootstrap_rows, sign_rows = ([], [], [])
    comparators = ['low_rank', *config['evaluation']['secondary_comparators']]
    for index, comparator in enumerate(comparators):
        result, bootstrap, signs = infer_pair(candidate, records[f'{comparator}_mse_top200'].to_numpy(float), comparator, int(config['evaluation']['primary_bootstrap_repetitions']), int(config['evaluation']['primary_bootstrap_seed']) + index, int(config['evaluation']['primary_sign_flip_repetitions']), int(config['evaluation']['primary_sign_flip_seed']) + index)
        comparisons.append(result)
        bootstrap_rows.extend(bootstrap)
        sign_rows.extend(signs)
    secondary_adjusted = transfer.holm_adjust([row['sign_flip_one_sided_p'] for row in comparisons if row['comparator'] != 'low_rank'])
    secondary_index = 0
    for row in comparisons:
        if row['comparator'] == 'low_rank':
            row['multiplicity_adjustment'] = 'single_predeclared_primary_none'
            row['adjusted_p'] = row['sign_flip_one_sided_p']
        else:
            row['multiplicity_adjustment'] = 'holm_across_secondary_comparators'
            row['adjusted_p'] = secondary_adjusted[secondary_index]
            secondary_index += 1
    comparison_map = {row['comparator']: row for row in comparisons}
    primary = comparison_map['low_rank']
    secondary_gene_sets = {name: float(1.0 - records[f'residual_program_mse_{name}'].mean() / records[f'low_rank_mse_{name}'].mean()) for name in ['top100', 'top500', 'all']}
    random_p = float((1 + np.count_nonzero(random_records['mean_mse_top200'] <= candidate.mean())) / (1 + len(random_records)))
    gate = config['promotion_gate']
    gates = {'minimum_improvement': primary['relative_improvement'] >= float(gate['minimum_relative_improvement_over_low_rank']), 'bootstrap_ci': primary['bootstrap_ci95'][0] > float(gate['primary_bootstrap_ci95_low_above']), 'sign_flip': primary['sign_flip_one_sided_p'] <= float(gate['primary_sign_flip_one_sided_p_at_most']), 'target_fraction': primary['targets_improved_fraction'] >= float(gate['targets_improved_fraction_at_least']), 'secondary_gene_sets': all((value > 0.0 for value in secondary_gene_sets.values())), 'random_program': random_p <= float(gate['random_program_exact_p_at_most'])}
    status = 'V30_LARGE_SCALE_RESIDUAL_PROGRAM_CONFIRMED' if all(gates.values()) else 'V30_LARGE_SCALE_RESIDUAL_PROGRAM_NOT_CONFIRMED'
    output_paths = {key: root_path(config['outputs'][key]) for key in ['record_metrics', 'random_program_metrics', 'bootstrap_metrics', 'sign_flip_metrics']}
    for path in output_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    records.to_csv(output_paths['record_metrics'], index=False)
    random_records.to_csv(output_paths['random_program_metrics'], index=False)
    pd.DataFrame(bootstrap_rows).to_csv(output_paths['bootstrap_metrics'], index=False)
    pd.DataFrame(sign_rows).to_csv(output_paths['sign_flip_metrics'], index=False)
    rankings = sorted([{'method': method, 'mean_mse_top200': float(records[f'{method}_mse_top200'].mean())} for method in ['residual_program', 'low_rank', 'gene_ridge', 'direct_program', 'scalar']], key=lambda row: row['mean_mse_top200'])
    summary = {'stage': 'V30 frozen large-scale internal development validation', 'generated_at_utc': datetime.now(timezone.utc).isoformat(), 'status': status, 'source_shape': list(map(int, source_data.shape)), 'target_shape': list(map(int, target_data.shape)), 'eligible_targets': len(eligible), 'common_measured_genes': len(genes), 'outer_fold_counts': records.groupby('outer_fold', observed=True).size().astype(int).to_dict(), 'primary_comparison': primary, 'comparisons': comparisons, 'secondary_gene_set_improvements_over_low_rank': secondary_gene_sets, 'random_program_exact_p': random_p, 'random_program_mean_mse_quantiles': np.quantile(random_records['mean_mse_top200'], [0.025, 0.5, 0.975]).tolist(), 'rankings': rankings, 'gates': {key: bool(value) for key, value in gates.items()}, 'claim_boundary': config['claim_boundary'], 'output_hashes': {key: digest(path) for key, path in output_paths.items()}}
    summary_path = root_path(config['outputs']['summary_json'])
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    root_path(config['outputs']['summary_markdown']).write_text('\n'.join(['# V30 large-scale residual-program validation', '', f'- Status: **{status}**', f'- Eligible targets: {len(eligible)}', f"- Residual program vs low-rank: {100 * primary['relative_improvement']:.3f}%", f"- Bootstrap 95% CI: [{100 * primary['bootstrap_ci95'][0]:.3f}%, {100 * primary['bootstrap_ci95'][1]:.3f}%]", f"- One-sided target sign-flip p: {primary['sign_flip_one_sided_p']:.6g}", f"- Targets improved: {100 * primary['targets_improved_fraction']:.2f}%", f"- Cohen dz: {primary['cohen_dz']:.4f}", f'- Random-program exact p: {random_p:.6f}', '', 'This is a same-study large-scale development validation, not a fresh external confirmation.', '']), encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False, indent=2))
if __name__ == '__main__':
    main()
