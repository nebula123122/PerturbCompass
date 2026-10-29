"""Check that registered protocols cannot silently substitute inputs."""
import json
from pathlib import Path

import pytest
import yaml

from perturbcompass import runtime


def test_required_input_missing_from_manifest_is_rejected(tmp_path, monkeypatch):
    (tmp_path / 'data').mkdir()
    (tmp_path / 'data/manifest.json').write_text(json.dumps({'files': []}), encoding='utf-8')
    monkeypatch.setattr(runtime, 'ROOT', tmp_path)
    with pytest.raises(RuntimeError, match='Unregistered raw input'):
        runtime.verify_raw_data(required={'data/raw/missing.h5ad'})


def test_duplicate_manifest_entries_are_rejected(tmp_path, monkeypatch):
    (tmp_path / 'data').mkdir()
    entry = {'path': 'data/raw/duplicate.h5ad'}
    (tmp_path / 'data/manifest.json').write_text(json.dumps({'files': [entry, entry]}), encoding='utf-8')
    monkeypatch.setattr(runtime, 'ROOT', tmp_path)
    with pytest.raises(RuntimeError, match='duplicate paths'):
        runtime.verify_raw_data()


def _fixture(tmp_path, monkeypatch):
    (tmp_path / 'docs/provenance').mkdir(parents=True)
    runner = tmp_path / 'analysis.py'; runner.write_text('pass\n', encoding='utf-8')
    config = tmp_path / 'protocol.yaml'; config.write_text('model: fixed\n', encoding='utf-8')
    files = {'analysis.py': runtime.digest(runner), 'protocol.yaml': runtime.digest(config)}
    (tmp_path / 'docs/provenance/reorganized_code_manifest.json').write_text(json.dumps({'files': files}), encoding='utf-8')
    monkeypatch.setattr(runtime, 'ROOT', tmp_path)
    monkeypatch.setattr(runtime, 'verify_small_artifacts', lambda: {})
    monkeypatch.setattr(runtime, 'verify_raw_data', lambda **kwargs: [])
    return runner, config


def test_external_configuration_is_rejected(tmp_path, monkeypatch):
    runner, registered = _fixture(tmp_path, monkeypatch)
    alternate = tmp_path / 'alternate.yaml'; alternate.write_text(registered.read_text(), encoding='utf-8')
    with pytest.raises(RuntimeError, match='Unregistered or modified analysis configuration'):
        runtime.verify_migrated_inputs(runner, alternate, {'model': 'fixed'})


def test_in_memory_protocol_change_is_rejected(tmp_path, monkeypatch):
    runner, config = _fixture(tmp_path, monkeypatch)
    with pytest.raises(RuntimeError, match='In-memory configuration'):
        runtime.verify_migrated_inputs(runner, config, {'model': 'changed'})


def test_registered_configuration_is_accepted(tmp_path, monkeypatch):
    runner, config = _fixture(tmp_path, monkeypatch)
    report = runtime.verify_migrated_inputs(runner, config, yaml.safe_load(config.read_text()))
    assert report['raw_inputs_verified'] == 0
