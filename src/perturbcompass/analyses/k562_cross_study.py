from __future__ import annotations
import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from perturbcompass.runtime import ROOT, root_path, digest, verify_migrated_inputs
from perturbcompass.transfer import parse_source_target, gene_lookup, clean_target_labels, dense_slice, aggregate_source, target_statistics, task_features, fit_weighted_ridge, scalar_oracles, program_source_features, fit_program_alphas, global_alpha, gene_wise_prediction, low_rank_prediction, mse, holm_adjust, infer_pair
from perturbcompass import transfer, programs, families, evaluation, apma, metadata
DEFAULT_CONFIG = ROOT / 'configs/k562_cross_study.yaml'

def verify_freeze(config_path, config):
    verify_migrated_inputs(__file__, config_path, config)
    return json.loads(root_path(config['outputs']['pre_run_manifest']).read_text(encoding='utf-8'))

def evaluate_source(source_name: str, source_data: ad.AnnData, target_data: dict[str, Any], union_gene_index: dict[str, int], config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    model_config = config['model']
    source_lookup = gene_lookup(source_data)
    genes = sorted(set(source_lookup) & set(union_gene_index))
    source_columns = np.asarray([source_lookup[gene] for gene in genes], dtype=np.int64)
    target_columns = np.asarray([union_gene_index[gene] for gene in genes], dtype=np.int64)
    source_targets, source_means, source_cells, source_constructs = aggregate_source(source_data, source_columns, int(config['source_preprocessing']['minimum_cells_per_construct']))
    source_map = {target: index for index, target in enumerate(source_targets)}
    eligible = sorted(set(source_targets) & set(target_data['target_means']) & set(genes))
    source = np.stack([source_means[source_map[target]] for target in eligible]).astype(np.float32)
    cells = np.asarray([source_cells[source_map[target]] for target in eligible], dtype=float)
    constructs = np.asarray([source_constructs[source_map[target]] for target in eligible], dtype=int)
    target_means = np.stack([target_data['target_means'][target][target_columns] for target in eligible]).astype(np.float32)
    target_cells = np.asarray([target_data['target_cells'][target] for target in eligible], dtype=int)
    minimum_targets = int(config['confirmation_gate']['minimum_eligible_targets'])
    if len(eligible) < minimum_targets:
        raise RuntimeError(f'{source_name} has only {len(eligible)} eligible targets; minimum is {minimum_targets}')
    program_count = int(model_config['program_count'])
    components = min(32, source.shape[0] - 1, source.shape[1] - 1)
    embedding = TruncatedSVD(n_components=components, random_state=int(model_config['program_seed'])).fit_transform(source.T)
    labels = KMeans(n_clusters=program_count, random_state=int(model_config['program_seed']), n_init=20).fit_predict(embedding)
    gene_to_index = {gene: index for index, gene in enumerate(genes)}
    target_programs = np.asarray([labels[gene_to_index[target]] for target in eligible], dtype=int)
    task_x = task_features(source, cells, constructs)
    signatures = [np.argsort(np.abs(source[index]))[-int(model_config['top_signature_genes']):] for index in range(len(eligible))]
    rng = np.random.default_rng(int(model_config['target_fold_seed']))
    folds = np.empty(len(eligible), dtype=int)
    permutation = rng.permutation(len(eligible))
    folds[permutation] = np.arange(len(eligible)) % int(model_config['outer_target_folds'])
    random_rng = np.random.default_rng(int(model_config['random_program_seed']))
    random_labels = [random_rng.permutation(labels) for _ in range(int(model_config['random_program_partitions']))]
    random_losses: dict[int, list[float]] = {index: [] for index in range(len(random_labels))}
    rows: list[dict[str, Any]] = []
    lower, upper = map(float, model_config['alpha_bounds'])
    for outer in range(int(model_config['outer_target_folds'])):
        train, test = (folds != outer, folds == outer)
        control_train = target_data['control_train_means'][outer][target_columns]
        control_test = target_data['control_test_means'][outer][target_columns]
        truth_train = target_means - control_train
        truth_test = target_means - control_test
        oracles, curvature = scalar_oracles(source, truth_train, signatures, lower, upper)
        scalar_prediction = np.clip(fit_weighted_ridge(task_x, oracles, curvature, train, test, float(model_config['scalar_ridge_alpha'])), lower, upper)
        program_alpha = fit_program_alphas(source, truth_train, train, test, labels, target_programs, task_x, scalar_prediction, model_config)
        fixed_alpha = global_alpha(source, truth_train, signatures, train, lower, upper)
        gene_prediction = gene_wise_prediction(source, truth_train, train, test, float(model_config['gene_wise_ridge_alpha']))
        low_rank = low_rank_prediction(source, truth_train, train, test, int(model_config['low_rank_components']), float(model_config['low_rank_ridge_alpha']), int(model_config['target_fold_seed']) + outer)
        test_indices = np.flatnonzero(test)
        random_alpha_by_partition = [fit_program_alphas(source, truth_train, train, test, random_partition, np.asarray([random_partition[gene_to_index[target]] for target in eligible], dtype=int), task_x, scalar_prediction, model_config) for random_partition in random_labels]
        for local_index, target_index in enumerate(test_indices):
            truth = truth_test[target_index]
            source_row = source[target_index]
            predictions = {'identity': source_row, 'global': fixed_alpha * source_row, 'scalar': scalar_prediction[local_index] * source_row, 'program': program_alpha[local_index][labels] * source_row, 'gene_ridge': gene_prediction[local_index], 'low_rank': low_rank[local_index]}
            ranked = np.argsort(np.abs(source_row))
            selections = {'top100': ranked[-min(100, len(ranked)):], 'top200': signatures[target_index], 'top500': ranked[-min(500, len(ranked)):], 'all': np.arange(len(genes))}
            row: dict[str, Any] = {'source': source_name, 'target': eligible[target_index], 'outer_fold': outer, 'source_cells': cells[target_index], 'source_constructs': constructs[target_index], 'target_cells': target_cells[target_index], 'scalar_alpha': scalar_prediction[local_index], 'program_alpha_mean': float(program_alpha[local_index].mean()), 'program_alpha_sd': float(program_alpha[local_index].std()), 'global_alpha': fixed_alpha}
            for selection_name, selected in selections.items():
                for method, prediction in predictions.items():
                    row[f'{method}_mse_{selection_name}'] = mse(prediction, truth, selected)
            rows.append(row)
            selected = selections['top200']
            for partition, alpha in enumerate(random_alpha_by_partition):
                prediction = alpha[local_index][random_labels[partition]] * source_row
                random_losses[partition].append(mse(prediction, truth, selected))
        print(json.dumps({'stage': 'v28_outer_fold', 'source': source_name, 'fold': outer, 'test_targets': int(test.sum())}), flush=True)
    records = pd.DataFrame(rows).sort_values('target').reset_index(drop=True)
    random_records = pd.DataFrame([{'source': source_name, 'partition': partition, 'mean_mse_top200': float(np.mean(losses))} for partition, losses in random_losses.items()])
    candidate_mean = float(records['program_mse_top200'].mean())
    random_p = float((1 + np.count_nonzero(random_records['mean_mse_top200'] <= candidate_mean)) / (1 + len(random_records)))
    metadata = {'eligible_targets': len(eligible), 'eligible_target_names': eligible, 'common_measured_genes': len(genes), 'control_fold_sizes': target_data['control_fold_sizes'], 'random_program_exact_p': random_p, 'random_program_mean_mse_quantiles': np.quantile(random_records['mean_mse_top200'], [0.025, 0.5, 0.975]).tolist()}
    return (records, random_records, metadata)

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=Path, default=DEFAULT_CONFIG)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text(encoding='utf-8'))
    verify_freeze(config_path, config)
    feasibility = json.loads(root_path(config['inputs']['metadata_feasibility']).read_text(encoding='utf-8'))
    if feasibility.get('target_expression_accessed') is not False:
        raise RuntimeError('metadata feasibility audit does not prove outcome blinding')
    target = ad.read_h5ad(root_path(config['inputs']['target_h5ad']))
    sources = {name: ad.read_h5ad(root_path(source['h5ad'])) for name, source in config['inputs']['sources'].items()}
    target_lookup = gene_lookup(target)
    if 'ncells' in target.var:
        retained_target_genes = {gene for gene, index in target_lookup.items() if float(target.var.iloc[index]['ncells']) >= int(config['target_preprocessing']['minimum_cells_expressing_gene'])}
    else:
        retained_target_genes = set(target_lookup)
    source_gene_union = set().union(*(set(gene_lookup(data)) for data in sources.values()))
    union_genes = sorted(retained_target_genes & source_gene_union)
    target_stats = target_statistics(target, union_genes, config['target_preprocessing'])
    union_gene_index = {gene: index for index, gene in enumerate(union_genes)}
    record_frames, random_frames = ([], [])
    source_metadata: dict[str, Any] = {}
    for source_name in [config['protocol']['primary_source'], config['protocol']['sensitivity_source']]:
        records, random_records, metadata = evaluate_source(source_name, sources[source_name], target_stats, union_gene_index, config)
        record_frames.append(records)
        random_frames.append(random_records)
        source_metadata[source_name] = metadata
    records = pd.concat(record_frames, ignore_index=True)
    random_records = pd.concat(random_frames, ignore_index=True)
    comparisons: list[dict[str, Any]] = []
    bootstrap_rows: list[dict[str, Any]] = []
    sign_rows: list[dict[str, Any]] = []
    comparators = ['identity', 'global', 'scalar', 'gene_ridge', 'low_rank']
    for source_index, source_name in enumerate([config['protocol']['primary_source'], config['protocol']['sensitivity_source']]):
        selected = records.loc[records['source'] == source_name].reset_index(drop=True)
        local = []
        for comparator_index, comparator in enumerate(comparators):
            result, boot, signs = infer_pair(selected['program_mse_top200'].to_numpy(dtype=float), selected[f'{comparator}_mse_top200'].to_numpy(dtype=float), source_name, comparator, config, 100 * source_index + comparator_index)
            local.append(result)
            bootstrap_rows.extend(boot)
            sign_rows.extend(signs)
        adjusted = holm_adjust([row['sign_flip_p'] for row in local])
        for row, value in zip(local, adjusted, strict=True):
            row['holm_adjusted_p'] = value
            comparisons.append(row)
    comparison_map = {(row['source'], row['comparator']): row for row in comparisons}
    primary = str(config['protocol']['primary_source'])
    sensitivity = str(config['protocol']['sensitivity_source'])
    primary_records = records.loc[records['source'] == primary]
    primary_scalar = comparison_map[primary, 'scalar']
    sensitivity_scalar = comparison_map[sensitivity, 'scalar']
    sensitivities = {name: float(1.0 - primary_records[f'program_mse_{name}'].mean() / primary_records[f'scalar_mse_{name}'].mean()) for name in ['top100', 'top500', 'all']}
    gate_config = config['confirmation_gate']
    gates = {'minimum_targets': source_metadata[primary]['eligible_targets'] >= int(gate_config['minimum_eligible_targets']), 'improvement': primary_scalar['improvement'] > float(gate_config['improvement_over_scalar_ridge_above']), 'bootstrap_ci': primary_scalar['ci95'][0] > float(gate_config['bootstrap_ci95_low_above']), 'holm_sign_flip': primary_scalar['holm_adjusted_p'] <= float(gate_config['holm_adjusted_sign_flip_p_at_most']), 'target_fraction': primary_scalar['targets_improved_fraction'] >= float(gate_config['targets_improved_fraction_at_least']), 'secondary_gene_sets': all((value > 0.0 for value in sensitivities.values())), 'random_program': source_metadata[primary]['random_program_exact_p'] <= float(gate_config['random_program_exact_p_at_most']), 'sensitivity_source': sensitivity_scalar['improvement'] > 0.0}
    publication_config = config['publication_superiority_gate']
    publication_gates = {'below_gene_ridge': comparison_map[primary, 'gene_ridge']['improvement'] > 0.0, 'below_low_rank': comparison_map[primary, 'low_rank']['improvement'] > 0.0, 'gene_ridge_holm': comparison_map[primary, 'gene_ridge']['holm_adjusted_p'] <= float(publication_config['holm_adjusted_p_at_most_for_each']), 'low_rank_holm': comparison_map[primary, 'low_rank']['holm_adjusted_p'] <= float(publication_config['holm_adjusted_p_at_most_for_each'])}
    core_status = 'CROSS_STUDY_CORE_CONFIRMED' if all(gates.values()) else 'CROSS_STUDY_CORE_NOT_CONFIRMED'
    publication_status = 'PUBLICATION_SUPERIORITY_CONFIRMED' if all(publication_gates.values()) else 'PUBLICATION_SUPERIORITY_NOT_CONFIRMED'
    output_paths = {key: root_path(config['outputs'][key]) for key in ['record_metrics', 'bootstrap_metrics', 'sign_flip_metrics', 'random_program_metrics']}
    for path in output_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    records.to_csv(output_paths['record_metrics'], index=False)
    pd.DataFrame(bootstrap_rows).to_csv(output_paths['bootstrap_metrics'], index=False)
    pd.DataFrame(sign_rows).to_csv(output_paths['sign_flip_metrics'], index=False)
    random_records.to_csv(output_paths['random_program_metrics'], index=False)
    rankings = {}
    for source_name in [primary, sensitivity]:
        selected = records.loc[records['source'] == source_name]
        rankings[source_name] = sorted([{'method': method, 'mean_mse_top200': float(selected[f'{method}_mse_top200'].mean())} for method in ['identity', 'global', 'scalar', 'program', 'gene_ridge', 'low_rank']], key=lambda row: row['mean_mse_top200'])
    summary = {'stage': 'V28 frozen Adamson cross-study K562 confirmation', 'generated_at_utc': datetime.now(timezone.utc).isoformat(), 'core_status': core_status, 'publication_status': publication_status, 'target_shape': list(map(int, target.shape)), 'target_qc_cells': target_stats['qc_cells'], 'target_control_cells': target_stats['control_cells'], 'source_metadata': source_metadata, 'comparisons': comparisons, 'primary_secondary_gene_set_improvements_over_scalar': sensitivities, 'rankings': rankings, 'gates': {key: bool(value) for key, value in gates.items()}, 'publication_gates': {key: bool(value) for key, value in publication_gates.items()}, 'claim_boundary': config['claim_boundary'], 'output_hashes': {key: digest(path) for key, path in output_paths.items()}}
    summary_path = root_path(config['outputs']['summary_json'])
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    primary_gene = comparison_map[primary, 'gene_ridge']
    primary_low_rank = comparison_map[primary, 'low_rank']
    root_path(config['outputs']['summary_markdown']).write_text('\n'.join(['# V28 Adamson cross-study K562 confirmation', '', f'- Core status: **{core_status}**', f'- Publication superiority: **{publication_status}**', f"- Primary eligible targets: {source_metadata[primary]['eligible_targets']}", f"- Program vs scalar improvement: {100 * primary_scalar['improvement']:.3f}%", f"- Bootstrap 95% CI: [{100 * primary_scalar['ci95'][0]:.3f}%, {100 * primary_scalar['ci95'][1]:.3f}%]", f"- Holm-adjusted sign-flip p: {primary_scalar['holm_adjusted_p']:.6f}", f"- Program vs gene-wise Ridge: {100 * primary_gene['improvement']:.3f}%", f"- Program vs low-rank Ridge: {100 * primary_low_rank['improvement']:.3f}%", f"- Exact-size random-program p: {source_metadata[primary]['random_program_exact_p']:.6f}", '', 'Training and test target effects use mutually exclusive target-control cells within every outer fold.', '']), encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False, indent=2))
if __name__ == '__main__':
    main()
