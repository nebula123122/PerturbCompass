"""Functional OOD and budget runs must verify the inputs they actually use."""
import json

import pytest
import yaml

from perturbcompass import runtime


@pytest.mark.parametrize('analysis', ['source_ood', 'anchor_budget'])
def test_sensitivities_verify_their_confirmation_source_and_target(analysis, monkeypatch):
    config_path = runtime.ROOT / 'configs' / (analysis + '.yaml')
    config = yaml.safe_load(config_path.read_text(encoding='utf-8'))
    observed = set()

    def capture(*, include_hashes=False, required=None):
        assert include_hashes is True
        observed.update(required)
        return [{'path': name} for name in required]

    monkeypatch.setattr(runtime, 'verify_raw_data', capture)
    runner = runtime.ROOT / 'src/perturbcompass/analyses' / (analysis + '.py')
    report = runtime.verify_migrated_inputs(runner, config_path, config)
    assert observed == {
        'data/raw/replogle_k562/K562_essential_normalized_bulk_01.h5ad',
        'data/raw/primary_cd4t/GWCD4i.DE_stats.h5ad',
    }
    assert report['raw_inputs_verified'] == 2
