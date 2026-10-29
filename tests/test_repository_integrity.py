"""Protect frozen evidence and the result identities during reorganization."""
from perturbcompass.reference_results import audit_reference_results
from perturbcompass.runtime import ROOT, verify_small_artifacts


def test_retained_evidence_hashes_match():
    counts = verify_small_artifacts()
    assert counts['immutable_artifacts.json'] > 0


def test_confirmation_and_full_anchor_budget_reconcile():
    report = audit_reference_results()
    assert report['status'] == 'REFERENCE_RESULTS_RECONCILED'


def test_public_file_names_use_scientific_purpose():
    import re
    for folder in ['src', 'tests', 'configs', 'docs', 'results/reference', 'data/processed']:
        for path in (ROOT / folder).rglob('*'):
            if path.is_file() and '__pycache__' not in path.parts:
                assert not re.search(r'(^|[/_.-])v\d+([/_.-]|$)', str(path.relative_to(ROOT)), re.I)
