"""Repository paths and post-result integrity checks.

Historical evidence hashes remain separate from the current code manifest.
Analysis entry points must execute registered code with their registered config.
"""
from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]


def root_path(value: str | Path) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def digest(path: Path, algorithm: str = 'sha256') -> str:
    value = hashlib.new(algorithm)
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(16 * 1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def verify_small_artifacts() -> dict[str, int]:
    counts = {}
    for name in ['reorganized_code_manifest.json', 'immutable_artifacts.json']:
        manifest = json.loads((ROOT / 'docs/provenance' / name).read_text(encoding='utf-8'))
        for relative, expected in manifest['files'].items():
            path = root_path(relative)
            if not path.is_relative_to(ROOT) or not path.is_file():
                raise RuntimeError(f'Missing repository artifact: {relative}')
            if digest(path) != expected:
                raise RuntimeError(f'Integrity mismatch: {relative}')
        counts[name] = len(manifest['files'])
    return counts


@lru_cache(maxsize=64)
def _verified_raw(relative: str, size: int, mtime_ns: int, expected: str) -> bool:
    if digest(root_path(relative)) != expected:
        raise RuntimeError(f'Raw data checksum mismatch: {relative}')
    return True


def verify_raw_data(*, include_hashes: bool = False, required: set[str] | None = None) -> list[dict[str, Any]]:
    manifest = json.loads((ROOT / 'data/manifest.json').read_text(encoding='utf-8'))
    entries = {entry['path']: entry for entry in manifest['files']}
    if len(entries) != len(manifest['files']):
        raise RuntimeError('Raw input manifest contains duplicate paths')
    if required is not None:
        unknown = required - set(entries)
        if unknown:
            raise RuntimeError('Unregistered raw input: ' + ', '.join(sorted(unknown)))
    report = []
    for relative, entry in entries.items():
        if required is not None and relative not in required:
            continue
        path = root_path(relative)
        if not path.is_relative_to(ROOT) or not path.is_file() or path.stat().st_size != entry['size_bytes']:
            raise RuntimeError(f'Missing or incomplete raw input: {relative}; see data/manifest.json')
        if include_hashes:
            stat = path.stat()
            _verified_raw(relative, stat.st_size, stat.st_mtime_ns, entry['checksums']['sha256'])
        report.append({'path': relative, 'size_verified': True, 'sha256_verified': include_hashes})
    return report


def _strings(value):
    if isinstance(value, dict):
        for item in value.values():
            yield from _strings(item)
    elif isinstance(value, list):
        for item in value:
            yield from _strings(item)
    elif isinstance(value, str):
        yield value


def _registered_file(value, files, role):
    path = root_path(value)
    if not path.is_relative_to(ROOT):
        raise RuntimeError(f'Unregistered {role}: {path}')
    relative = path.relative_to(ROOT).as_posix()
    if relative not in files or not path.is_file() or digest(path) != files[relative]:
        raise RuntimeError(f'Unregistered or modified {role}: {relative}')
    return path


def verify_migrated_inputs(runner, config_path, config):
    import yaml

    verify_small_artifacts()
    files = json.loads((ROOT / 'docs/provenance/reorganized_code_manifest.json').read_text(encoding='utf-8'))['files']
    _registered_file(runner, files, 'analysis runner')
    config_file = _registered_file(config_path, files, 'analysis configuration')
    if yaml.safe_load(config_file.read_text(encoding='utf-8')) != config:
        raise RuntimeError('In-memory configuration differs from its registered file')
    required = set()
    seen = set()

    def inspect(value):
        for token in _strings(value):
            if token.startswith('data/raw/'):
                required.add(token)
            elif token.startswith('configs/') and token.endswith('.yaml') and token not in seen:
                seen.add(token)
                nested = _registered_file(token, files, 'parent configuration')
                inspect(yaml.safe_load(nested.read_text(encoding='utf-8')))

    inspect(config)
    checked = verify_raw_data(include_hashes=True, required=required)
    return {'classification': 'post_result_refactor_integrity_check', 'raw_inputs_verified': len(checked)}
