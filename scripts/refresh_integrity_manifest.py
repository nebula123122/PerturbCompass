"""Refresh current code hashes after reviewed edits; never rewrite frozen evidence."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
DIRECTORIES = ['src', 'tests', 'configs', 'scripts', '.github']
ROOT_FILES = ['pyproject.toml', 'environment.yml', '.gitattributes']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reason', required=True, help='Explain the reviewed code or configuration change')
    args = parser.parse_args()
    if not args.reason.strip():
        parser.error('A nonempty change reason is required')
    manifest_path = ROOT / 'docs/provenance/reorganized_code_manifest.json'
    old = json.loads(manifest_path.read_text(encoding='utf-8'))
    files = {}
    for directory in DIRECTORIES:
        for path in (ROOT / directory).rglob('*'):
            if path.is_file() and not any(part == '__pycache__' or part.endswith('.egg-info') for part in path.parts):
                files[path.relative_to(ROOT).as_posix()] = digest(path)
    for name in ROOT_FILES:
        if (ROOT / name).is_file():
            files[name] = digest(ROOT / name)
    before = old['files']
    changed = {name: {'before': before.get(name), 'after': files.get(name)} for name in sorted(set(before) | set(files)) if before.get(name) != files.get(name)}
    if not changed:
        print('No registered code changes.'); return
    timestamp = datetime.now(timezone.utc).isoformat()
    history_path = ROOT / 'docs/provenance/code_revisions.json'
    history = json.loads(history_path.read_text(encoding='utf-8')) if history_path.exists() else {'classification': 'post_result_code_maintenance_not_preregistration', 'revisions': []}
    head = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, capture_output=True, text=True)
    history['revisions'].append({'recorded_at_utc': timestamp, 'reason': args.reason, 'base_commit': head.stdout.strip() if head.returncode == 0 else None, 'previous_manifest_sha256': digest(manifest_path), 'changed_files': changed, 'historical_evidence_modified': False})
    history_path.write_text(json.dumps(history, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    manifest_path.write_text(json.dumps({'generated_at_utc': timestamp, 'classification': 'post_result_code_integrity_manifest', 'files': dict(sorted(files.items()))}, indent=2) + '\n', encoding='utf-8')
    print(f'Registered {len(files)} files; recorded {len(changed)} code/configuration changes.')


if __name__ == '__main__':
    main()
