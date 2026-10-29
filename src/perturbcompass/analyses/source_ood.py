from __future__ import annotations
from perturbcompass import families as family_models
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
import anndata as ad
import h5py
import numpy as np
import pandas as pd
import yaml
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD
from perturbcompass.runtime import ROOT, root_path, digest, verify_migrated_inputs
from perturbcompass import transfer, programs, families, evaluation, apma, metadata
CONFIG = ROOT / 'configs/source_ood.yaml'
CONFIRMATION_CONFIG = ROOT / 'configs/cd4t_confirmation.yaml'
SEAL = ROOT / 'docs/provenance/source_ood/source_ood_pre_run_manifest.json'

def verify_seal(design):
    verify_migrated_inputs(__file__, CONFIG, design)
    return json.loads(root_path(design['outputs']['pre_numeric_manifest']).read_text(encoding='utf-8'))

def stratified_bootstrap_regret(selector_loss: np.ndarray, family_losses: np.ndarray, clusters: np.ndarray, repetitions: int, seed: int) -> tuple[np.ndarray, float, list[float]]:
    actual = float(selector_loss.mean() / family_losses.mean(axis=0).min() - 1.0)
    groups = [np.flatnonzero(clusters == value) for value in sorted(np.unique(clusters))]
    rng = np.random.default_rng(seed)
    values = np.empty(repetitions, dtype=float)
    for repetition in range(repetitions):
        rows = np.concatenate([rng.choice(group, size=len(group), replace=True) for group in groups])
        values[repetition] = selector_loss[rows].mean() / family_losses[rows].mean(axis=0).min() - 1.0
    interval = np.quantile(values, [0.025, 0.975])
    return (values, actual, [float(interval[0]), float(interval[1])])

def main() -> None:
    design = yaml.safe_load(CONFIG.read_text(encoding='utf-8'))
    confirmation = yaml.safe_load(CONFIRMATION_CONFIG.read_text(encoding='utf-8'))
    verify_seal(design)
    partition = pd.read_csv(root_path(design['outputs']['source_cluster_partition'])).sort_values('target').reset_index(drop=True)
    partition['target'] = partition['target'].astype(str).str.upper()
    targets = partition['target'].tolist()
    anchors_all = partition['role'].eq('anchor').to_numpy()
    tests_all = partition['role'].eq('test').to_numpy()
    cluster_labels = partition['source_ood_cluster'].to_numpy(dtype=int)
    target_path = root_path(confirmation['target']['path'])
    with h5py.File(target_path, 'r') as handle:
        var = handle['var']
        measured = metadata.read_frame_column(var, 'gene_name') if 'gene_name' in var else metadata.read_index(var)
    target_genes = {str(value).strip().upper() for value in measured if str(value).strip()}
    source_data = ad.read_h5ad(root_path(confirmation['source']['h5ad']))
    source_lookup = transfer.gene_lookup(source_data)
    genes = sorted(set(source_lookup) & target_genes)
    source_columns = np.asarray([source_lookup[gene] for gene in genes], dtype=np.int64)
    source_targets, source_means, source_cells, source_constructs = transfer.aggregate_source(source_data, source_columns, int(confirmation['source']['minimum_cells_per_construct']))
    del source_data
    source_map = {str(target).upper(): index for index, target in enumerate(source_targets)}
    source = np.stack([source_means[source_map[target]] for target in targets]).astype(np.float32)
    cells = np.asarray([source_cells[source_map[target]] for target in targets], dtype=float)
    constructs = np.asarray([source_constructs[source_map[target]] for target in targets], dtype=int)
    source_rms = np.sqrt(np.maximum(np.mean(source * source, axis=1), 1e-12))
    model = confirmation['model']
    embedding = TruncatedSVD(n_components=min(32, len(targets) - 1, len(genes) - 1), random_state=int(model['program_seed'])).fit_transform(source.T)
    gene_programs = KMeans(n_clusters=int(model['program_count']), random_state=int(model['program_seed']), n_init=20).fit_predict(embedding)
    target_programs, fallback_targets = apma.assign_target_programs(targets, genes, gene_programs, source)
    task_x = transfer.task_features(source, cells, constructs)
    signatures = evaluation.scope_indices(source, int(model['top_signature_genes']))
    families = list(map(str, model['families']))
    anchor_rows: list[dict[str, Any]] = []
    heldout_rows: list[dict[str, Any]] = []
    bootstrap_rows: list[dict[str, Any]] = []
    summaries: dict[str, Any] = {}
    for condition_offset, condition in enumerate(confirmation['target']['conditions']):
        truth, _ = apma.target_log_fc(target_path, genes, targets, str(condition), str(confirmation['target']['effect_layer']))
        condition_rows = []
        cluster_summaries = {}
        for cluster in sorted(np.unique(cluster_labels)):
            local_anchors = anchors_all & (cluster_labels != cluster)
            local_test = tests_all & (cluster_labels == cluster)
            anchor_indices = np.flatnonzero(local_anchors)
            folds = programs.stratified_round_robin_folds(source_rms[anchor_indices], int(model['anchor_oof_folds']), int(model['anchor_oof_seed']) + int(cluster))
            anchor_loss = np.full((len(anchor_indices), len(families)), np.nan)
            positions = {int(index): pos for pos, index in enumerate(anchor_indices)}
            for fold in range(int(model['anchor_oof_folds'])):
                validation_indices = anchor_indices[folds == fold]
                validation = np.zeros(len(targets), dtype=bool)
                validation[validation_indices] = True
                train = local_anchors & ~validation
                predictions = family_models.family_predictions(source, truth, train, validation, source_rms, gene_programs, target_programs, task_x, signatures, model, 3000 + 100 * int(cluster) + fold)
                for local_position, target_index in enumerate(validation_indices):
                    selected = signatures[target_index]
                    row: dict[str, Any] = {'condition': condition, 'cluster': int(cluster), 'target': targets[target_index], 'fold': int(fold)}
                    for family_index, family in enumerate(families):
                        loss = transfer.mse(predictions[family][local_position], truth[target_index], selected)
                        anchor_loss[positions[int(target_index)], family_index] = loss
                        row[f'{family}_mse_top200'] = loss
                    anchor_rows.append(row)
            decision = apma.choose_probability_selector(anchor_loss, families, int(design['selector']['anchor_bootstrap_repetitions']), 380010 + 100 * condition_offset + int(cluster), float(design['selector']['dominant_probability_threshold']))
            predictions = family_models.family_predictions(source, truth, local_anchors, local_test, source_rms, gene_programs, target_programs, task_x, signatures, model, 4000 + 100 * condition_offset + int(cluster))
            selected_prediction = apma.weighted_prediction(predictions, families, np.asarray(decision['weight_vector']))
            test_indices = np.flatnonzero(local_test)
            local_signatures = [signatures[index] for index in test_indices]
            truth_test = truth[local_test]
            family_loss = {family: evaluation.per_target_mse(predictions[family], truth_test, local_signatures) for family in families}
            selector_loss = evaluation.per_target_mse(selected_prediction, truth_test, local_signatures)
            residual_improvement = float(1.0 - family_loss['residual_program'].mean() / family_loss['low_rank'].mean())
            cluster_summaries[str(cluster)] = {'test_targets': int(local_test.sum()), 'training_anchors': int(local_anchors.sum()), 'selector': {key: value for key, value in decision.items() if key != 'weight_vector'}, 'residual_vs_low_rank_relative_improvement': residual_improvement, 'residual_beats_low_rank': bool(residual_improvement > 0)}
            for local_position, target_index in enumerate(test_indices):
                row = {'condition': condition, 'cluster': int(cluster), 'target': targets[target_index], 'selector_mse_top200': float(selector_loss[local_position])}
                for family in families:
                    row[f'{family}_mse_top200'] = float(family_loss[family][local_position])
                condition_rows.append(row)
                heldout_rows.append(row)
        frame = pd.DataFrame(condition_rows).sort_values(['cluster', 'target']).reset_index(drop=True)
        selector_loss = frame['selector_mse_top200'].to_numpy(dtype=float)
        family_matrix = np.column_stack([frame[f'{family}_mse_top200'].to_numpy(dtype=float) for family in families])
        bootstrap, regret, interval = stratified_bootstrap_regret(selector_loss, family_matrix, frame['cluster'].to_numpy(dtype=int), int(design['evaluation']['bootstrap_repetitions']), int(design['evaluation']['bootstrap_seed']) + condition_offset)
        for repetition, value in enumerate(bootstrap):
            bootstrap_rows.append({'condition': condition, 'bootstrap': repetition, 'regret': value})
        positive_clusters = sum((value['residual_beats_low_rank'] for value in cluster_summaries.values()))
        summaries[str(condition)] = {'clusters': cluster_summaries, 'heldout_targets': len(frame), 'selector_regret': regret, 'selector_regret_cluster_stratified_bootstrap_ci95': interval, 'regret_gate_upper_ci_at_most_0_05': bool(interval[1] <= 0.05), 'residual_positive_clusters': int(positive_clusters), 'residual_rule_at_least_six_clusters': bool(positive_clusters >= 6), 'fixed_family_mean_mse_top200': {family: float(frame[f'{family}_mse_top200'].mean()) for family in families}}
        del truth
    primary_pass = all((value['regret_gate_upper_ci_at_most_0_05'] for value in summaries.values()))
    status = 'V38_SOURCE_OOD_SENSITIVITY_SUPPORTS_SELECTOR' if primary_pass else 'V38_SOURCE_OOD_SENSITIVITY_DOES_NOT_SUPPORT_SELECTOR'
    summary = {'status': status, 'generated_at_utc': datetime.now(timezone.utc).isoformat(), 'conditions': summaries, 'actual_source_cluster_count': int(len(np.unique(cluster_labels))), 'cluster_feasibility_note': 'Original 8-cluster layout failed the frozen 40-test-target minimum; the largest source-only feasible candidate was used before target numeric access.', 'joint_stress_test_gate': {'passed': primary_pass, 'rule': design['evaluation']['stress_test_support_rule']}, 'target_program_assignment': {'method': 'measured-target gene program when available; otherwise the program with highest mean squared source-effect energy', 'fallback_target_count': len(fallback_targets), 'fallback_targets': fallback_targets}, 'claim_boundary': design['claim_boundary']}
    output_dir = ROOT / 'results/reproduced/source_ood'
    output_dir.mkdir(parents=True, exist_ok=True)
    anchor_path = output_dir / 'anchor_oof_metrics.csv'
    heldout_path = output_dir / 'heldout_metrics.csv'
    bootstrap_path = output_dir / 'bootstrap_regret.csv'
    pd.DataFrame(anchor_rows).to_csv(anchor_path, index=False)
    pd.DataFrame(heldout_rows).to_csv(heldout_path, index=False)
    pd.DataFrame(bootstrap_rows).to_csv(bootstrap_path, index=False)
    summary_path = root_path(design['outputs']['summary_json'])
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    lines = ['# V38 source-defined OOD 压力测试', '', f'- 状态：**{status}**', f"- 三状态 5% regret 联合门槛：**{('通过' if primary_pass else '未通过')}**。", '']
    for condition, value in summaries.items():
        lines += [f'## {condition}', '', f"- regret：{100 * value['selector_regret']:.3f}%（95% CI {100 * value['selector_regret_cluster_stratified_bootstrap_ci95'][0]:.3f}% 至 {100 * value['selector_regret_cluster_stratified_bootstrap_ci95'][1]:.3f}%）", f"- residual 胜过 low-rank 的 source clusters：{value['residual_positive_clusters']}/{len(np.unique(cluster_labels))}", '']
    lines.append('该分析只支持 source geometry 外推敏感性，不等同于按目标真值聚类的 EL-PERTURB OOD，也不构成新的独立研究。')
    root_path(design['outputs']['summary_markdown']).write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False, indent=2))
if __name__ == '__main__':
    main()
