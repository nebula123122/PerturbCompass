from __future__ import annotations
from perturbcompass.analyses import feng_external_validation as feng_data
from perturbcompass import transfer, programs, evaluation, apma, design as prepare
import gzip
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

def write_lfc(path: Path, rows: list[dict[str, object]]) -> None:
    frame = pd.DataFrame(rows)
    with gzip.open(path, 'wt', encoding='utf-8', newline='') as handle:
        frame.to_csv(handle, sep='\t', index=False)

def test_load_lfc_matrix_preserves_frozen_target_order(tmp_path: Path) -> None:
    path = tmp_path / 'lfc.tsv.gz'
    rows = []
    for target, offset in [('B', 2.0), ('A', 1.0)]:
        for gene, value in [('g2', 0.2), ('g1', 0.1)]:
            rows.append({'Target': target, 'Expressed_Gene_Symbol': gene, 'Expressed_Gene_Ens_ID': f'ENSG_{gene}', 'lfc': value + offset, 'pval_adj': 0.5})
    write_lfc(path, rows)
    matrix = feng_data.load_lfc_matrix(path, ['A', 'B'], chunk_size=2)
    assert matrix.index.tolist() == ['A', 'B']
    assert matrix.columns.tolist() == ['G1', 'G2']
    assert np.allclose(matrix.to_numpy(), [[1.1, 1.2], [2.1, 2.2]])

def test_load_lfc_matrix_rejects_duplicate_target_gene(tmp_path: Path) -> None:
    path = tmp_path / 'duplicate.tsv.gz'
    row = {'Target': 'A', 'Expressed_Gene_Symbol': 'G1', 'Expressed_Gene_Ens_ID': 'ENSG1', 'lfc': 0.1, 'pval_adj': 0.5}
    write_lfc(path, [row, row])
    with pytest.raises(RuntimeError, match='duplicate target-gene'):
        feng_data.load_lfc_matrix(path, ['A'], chunk_size=1)

def test_status_requires_replication_and_publication_guardrail() -> None:
    replication = {'a': True, 'b': True}
    publication = {'replication': True, 'ridge': True}
    assert feng_data.status_from_gates(replication, publication) == 'V31_EXTERNAL_REPLICATION_AND_GENE_RIDGE_NONINFERIOR'
    publication['ridge'] = False
    assert feng_data.status_from_gates(replication, publication) == 'V31_EXTERNAL_LOW_RANK_REPLICATION_ONLY'
    replication['a'] = False
    assert feng_data.status_from_gates(replication, publication) == 'V31_EXTERNAL_REPLICATION_NOT_CONFIRMED'
