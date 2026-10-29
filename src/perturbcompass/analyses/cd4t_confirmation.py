from __future__ import annotations
from perturbcompass import families as family_models
import argparse
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
from perturbcompass.apma import choose_probability_selector, weighted_prediction, assign_target_programs, target_log_fc
from perturbcompass import transfer, programs, families, evaluation, apma, metadata
DEFAULT_CONFIG = ROOT / 'configs/cd4t_confirmation.yaml'
COMPATIBILITY_AMENDMENT = ROOT / 'docs/provenance/cd4t_confirmation/target_program_compatibility_amendment.json'

def verify_seal(config_path, config):
    verify_migrated_inputs(__file__, config_path, config)
    return json.loads(root_path(config['outputs']['pre_numeric_manifest']).read_text(encoding='utf-8'))

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=Path, default=DEFAULT_CONFIG)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text(encoding='utf-8'))
    manifest = verify_seal(config_path, config)
    partition = pd.read_csv(root_path(config['frozen_inputs']['shared_partition'])).sort_values('target').reset_index(drop=True)
    partition['target'] = partition['target'].astype(str).str.upper()
    targets = partition['target'].tolist()
    anchors = partition['role'].eq(config['partition']['anchor_role']).to_numpy()
    heldout = partition['role'].eq(config['partition']['heldout_role']).to_numpy()
    if int(anchors.sum()) != int(config['partition']['expected_anchor_targets']):
        raise RuntimeError('anchor count changed')
    if int(heldout.sum()) != int(config['partition']['expected_heldout_targets']):
        raise RuntimeError('heldout count changed')
    target_path = root_path(config['target']['path'])
    with h5py.File(target_path, 'r') as handle:
        var = handle['var']
        measured = metadata.read_frame_column(var, 'gene_name') if 'gene_name' in var else metadata.read_index(var)
    target_genes = {str(value).strip().upper() for value in measured if str(value).strip()}
    source_data = ad.read_h5ad(root_path(config['source']['h5ad']))
    source_lookup = transfer.gene_lookup(source_data)
    genes = sorted(set(source_lookup) & target_genes)
    if len(genes) != int(config['target']['expected_common_genes']):
        raise RuntimeError(f"expected {config['target']['expected_common_genes']} common genes, found {len(genes)}")
    source_columns = np.asarray([source_lookup[gene] for gene in genes], dtype=np.int64)
    source_targets, source_means, source_cells, source_constructs = transfer.aggregate_source(source_data, source_columns, int(config['source']['minimum_cells_per_construct']))
    del source_data
    source_map = {str(target).upper(): i for i, target in enumerate(source_targets)}
    if set(targets) - set(source_map):
        raise RuntimeError('frozen targets disappeared from source')
    source = np.stack([source_means[source_map[target]] for target in targets]).astype(np.float32)
    cells = np.asarray([source_cells[source_map[target]] for target in targets], dtype=float)
    constructs = np.asarray([source_constructs[source_map[target]] for target in targets], dtype=int)
    source_rms = np.sqrt(np.maximum(np.mean(source * source, axis=1), 1e-12))
    model = config['model']
    components = min(32, len(targets) - 1, len(genes) - 1)
    embedding = TruncatedSVD(n_components=components, random_state=int(model['program_seed'])).fit_transform(source.T)
    labels = KMeans(n_clusters=int(model['program_count']), random_state=int(model['program_seed']), n_init=20).fit_predict(embedding)
    target_programs, fallback_targets = assign_target_programs(targets, genes, labels, source)
    task_x = transfer.task_features(source, cells, constructs)
    signatures = evaluation.scope_indices(source, int(model['top_signature_genes']))
    scopes = {'top100': evaluation.scope_indices(source, 100), 'top200': signatures, 'top500': evaluation.scope_indices(source, 500), 'all': evaluation.scope_indices(source, None)}
    families = list(map(str, model['families']))
    anchor_indices = np.flatnonzero(anchors)
    folds = programs.stratified_round_robin_folds(source_rms[anchor_indices], int(model['anchor_oof_folds']), int(model['anchor_oof_seed']))
    anchor_rows: list[dict[str, Any]] = []
    heldout_rows: list[dict[str, Any]] = []
    bootstrap_rows: list[dict[str, Any]] = []
    summaries: dict[str, Any] = {}
    comparisons: dict[str, float] = {}
    for offset, condition in enumerate(config['target']['conditions']):
        truth, data_report = target_log_fc(target_path, genes, targets, str(condition), str(config['target']['effect_layer']))
        anchor_loss = np.full((len(anchor_indices), len(families)), np.nan)
        positions = {int(index): pos for pos, index in enumerate(anchor_indices)}
        for fold in range(int(model['anchor_oof_folds'])):
            test_indices = anchor_indices[folds == fold]
            local_test = np.zeros(len(targets), dtype=bool)
            local_test[test_indices] = True
            local_train = anchors & ~local_test
            predictions = family_models.family_predictions(source, truth, local_train, local_test, source_rms, labels, target_programs, task_x, signatures, model, fold)
            for local_pos, target_index in enumerate(test_indices):
                selected = signatures[target_index]
                row: dict[str, Any] = {'condition': condition, 'target': targets[target_index], 'fold': int(fold)}
                for family_index, family in enumerate(families):
                    loss = transfer.mse(predictions[family][local_pos], truth[target_index], selected)
                    anchor_loss[positions[int(target_index)], family_index] = loss
                    row[f'{family}_mse_top200'] = loss
                anchor_rows.append(row)
        decision = choose_probability_selector(anchor_loss, families, int(config['selector']['anchor_bootstrap_repetitions']), int(config['selector']['anchor_bootstrap_seed']) + offset, float(config['selector']['dominant_probability_threshold']))
        predictions = family_models.family_predictions(source, truth, anchors, heldout, source_rms, labels, target_programs, task_x, signatures, model, 999 + offset)
        selector_prediction = weighted_prediction(predictions, families, np.asarray(decision['weight_vector']))
        heldout_indices = np.flatnonzero(heldout)
        truth_test = truth[heldout]
        local_scopes = {name: [values[index] for index in heldout_indices] for name, values in scopes.items()}
        loss_by_scope: dict[str, dict[str, np.ndarray]] = {}
        for scope, indices in local_scopes.items():
            loss_by_scope[scope] = {family: evaluation.per_target_mse(predictions[family], truth_test, indices) for family in families}
            loss_by_scope[scope]['selector'] = evaluation.per_target_mse(selector_prediction, truth_test, indices)
        primary_matrix = np.column_stack([loss_by_scope['top200'][family] for family in families])
        selector_loss = loss_by_scope['top200']['selector']
        regret_values, regret, regret_ci = evaluation.bootstrap_regret(selector_loss, primary_matrix, int(config['evaluation']['bootstrap_repetitions']), int(config['evaluation']['bootstrap_seed']) + offset)
        for repeat, value in enumerate(regret_values):
            bootstrap_rows.append({'condition': condition, 'bootstrap': repeat, 'regret': value})
        for local_index, target_index in enumerate(heldout_indices):
            row = {'condition': condition, 'target': targets[target_index]}
            for scope, losses in loss_by_scope.items():
                for family, values in losses.items():
                    row[f'{family}_mse_{scope}'] = float(values[local_index])
            heldout_rows.append(row)
        scalar_test = evaluation.paired_inference(selector_loss, loss_by_scope['top200']['scalar'], int(config['evaluation']['bootstrap_repetitions']), int(config['evaluation']['bootstrap_seed']) + 100 + offset, int(config['evaluation']['sign_flip_repetitions']), int(config['evaluation']['sign_flip_seed']) + offset)
        residual_test = evaluation.paired_inference(loss_by_scope['top200']['residual_program'], loss_by_scope['top200']['low_rank'], int(config['evaluation']['bootstrap_repetitions']), int(config['evaluation']['bootstrap_seed']) + 200 + offset, int(config['evaluation']['sign_flip_repetitions']), int(config['evaluation']['sign_flip_seed']) + 100 + offset)
        comparisons[f'{condition}:selector_vs_scalar'] = scalar_test['sign_flip_one_sided_p']
        comparisons[f'{condition}:residual_vs_low_rank'] = residual_test['sign_flip_one_sided_p']
        means = {family: float(loss_by_scope['top200'][family].mean()) for family in families}
        summaries[str(condition)] = {'data': data_report, 'selector': {key: value for key, value in decision.items() if key != 'weight_vector'}, 'heldout': {'targets': int(heldout.sum()), 'selector_mean_mse_top200': float(selector_loss.mean()), 'fixed_family_mean_mse_top200': means, 'oracle_fixed_family': min(means, key=means.get), 'selector_regret': regret, 'selector_regret_bootstrap_ci95': regret_ci, 'primary_gate_upper_ci_at_most_0_02': bool(regret_ci[1] <= float(config['evaluation']['maximum_regret_ci_upper'])), 'selector_vs_scalar': scalar_test, 'residual_program_vs_low_rank': residual_test, 'scope_relative_improvement_vs_scalar': {scope: float(1.0 - losses['selector'].mean() / losses['scalar'].mean()) for scope, losses in loss_by_scope.items()}}}
        del truth, predictions, selector_prediction
    adjusted = evaluation.holm_adjust(comparisons)
    for key, value in adjusted.items():
        condition, comparison = key.split(':')
        field = 'selector_vs_scalar' if comparison == 'selector_vs_scalar' else 'residual_program_vs_low_rank'
        summaries[condition]['heldout'][field]['holm_adjusted_p'] = float(value)
    primary_pass = all((summaries[str(condition)]['heldout']['primary_gate_upper_ci_at_most_0_02'] for condition in config['evaluation']['required_conditions']))
    status = 'V36_PROSPECTIVE_SELECTOR_CONFIRMED_ALL_THREE_CONTEXTS' if primary_pass else 'V36_PROSPECTIVE_SELECTOR_NOT_CONFIRMED_ALL_THREE_CONTEXTS'
    summary = {'status': status, 'generated_at_utc': datetime.now(timezone.utc).isoformat(), 'pre_numeric_manifest_sha256': digest(root_path(config['outputs']['pre_numeric_manifest'])), 'common_genes': len(genes), 'anchor_targets': int(anchors.sum()), 'heldout_targets': int(heldout.sum()), 'target_program_assignment': {'direct_gene_membership': len(targets) - len(fallback_targets), 'source_effect_energy_fallback': len(fallback_targets), 'fallback_targets': fallback_targets, 'rule': 'missing target genes use the program with maximum mean squared source-effect energy'}, 'conditions': summaries, 'primary_joint_gate': {'rule': config['evaluation']['primary_success_rule'], 'passed': primary_pass}, 'claim_boundary': config['claim_boundary']}
    paths = {key: root_path(value) for key, value in config['outputs'].items()}
    for key in ('anchor_oof_metrics', 'heldout_metrics', 'bootstrap_regret', 'summary_json', 'summary_markdown', 'output_audit'):
        paths[key].parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(anchor_rows).to_csv(paths['anchor_oof_metrics'], index=False)
    pd.DataFrame(heldout_rows).to_csv(paths['heldout_metrics'], index=False)
    pd.DataFrame(bootstrap_rows).to_csv(paths['bootstrap_regret'], index=False)
    paths['summary_json'].write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    lines = ['# V36 原代 CD4 T 细胞三状态前瞻确认', '', f'- 状态：**{status}**', f"- 联合门槛：**{('通过' if primary_pass else '未通过')}**（三种状态 regret 95% CI 上界均 ≤2%）。", '']
    for condition, value in summaries.items():
        held = value['heldout']
        selector = value['selector']
        weights = ', '.join((f'{name}={weight:.3f}' for name, weight in selector['weights'].items()))
        lines += [f'## {condition}', '', f"- anchor 决策：{selector['mode']}；权重：{weights}", f"- anchor 赢家：{selector['anchor_winner']}；稳定率：{100 * selector['leading_probability']:.2f}%", f"- held-out 最佳固定家族：{held['oracle_fixed_family']}", f"- selector regret：{100 * held['selector_regret']:.3f}%（95% CI {100 * held['selector_regret_bootstrap_ci95'][0]:.3f}% 至 {100 * held['selector_regret_bootstrap_ci95'][1]:.3f}%）", f"- selector 相对 scalar：{100 * held['selector_vs_scalar']['relative_improvement']:.3f}%", f"- residual 相对 low-rank：{100 * held['residual_program_vs_low_rank']['relative_improvement']:.3f}%", '']
    lines.append('三种状态和四位供者来自同一项研究，因此联合门槛检验跨状态稳定性，不等于三项独立外部复制。')
    paths['summary_markdown'].write_text('\n'.join(lines) + '\n', encoding='utf-8')
    output_files = [paths['anchor_oof_metrics'], paths['heldout_metrics'], paths['bootstrap_regret'], paths['summary_json'], paths['summary_markdown']]
    audit = {'status': 'V36_OUTPUT_AUDIT_PASS', 'summary_status': status, 'row_counts': {'anchor_oof': len(anchor_rows), 'heldout': len(heldout_rows), 'bootstrap_regret': len(bootstrap_rows)}, 'expected_row_counts': {'anchor_oof': 3 * int(anchors.sum()), 'heldout': 3 * int(heldout.sum()), 'bootstrap_regret': 3 * int(config['evaluation']['bootstrap_repetitions'])}, 'files': {str(path.relative_to(ROOT)).replace('\\', '/'): digest(path) for path in output_files}}
    if audit['row_counts'] != audit['expected_row_counts']:
        raise RuntimeError('V36 output row count audit failed')
    paths['output_audit'].write_text(json.dumps(audit, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False, indent=2))
if __name__ == '__main__':
    main()
