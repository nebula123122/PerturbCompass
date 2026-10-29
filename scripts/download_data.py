"""Download public paper inputs and verify size and SHA-256 before acceptance."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parents[1]


def checksum(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(16 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--list', action='store_true')
    parser.add_argument('--dataset', help='Study directory name from the manifest')
    args = parser.parse_args()
    entries = json.loads((ROOT / 'data/manifest.json').read_text())['files']
    if args.list:
        for entry in entries:
            print(entry['path'], entry['size_bytes'], entry.get('download_url', 'manual author artifact'))
        return
    if not args.dataset:
        parser.error('Use --list or --dataset STUDY_NAME')
    selected = [entry for entry in entries if Path(entry['path']).parts[2] == args.dataset]
    if not selected:
        parser.error('Unknown study directory; use --list')
    for entry in selected:
        destination = ROOT / entry['path']
        if destination.exists():
            if destination.stat().st_size == entry['size_bytes'] and checksum(destination) == entry['checksums']['sha256']:
                print('Verified existing', entry['path']); continue
            raise RuntimeError(f'Existing input differs: {destination}; inspect it before retrying.')
        if 'download_url' not in entry:
            raise RuntimeError(f'No byte-verified public endpoint recorded for {entry["path"]}; obtain the author artifact specified in data/manifest.json.')
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_name(destination.name + '.download')
        if temporary.exists():
            raise RuntimeError(f'Partial download exists: {temporary}; inspect before retrying.')
        print('Downloading', entry['path'], flush=True)
        request = urllib.request.Request(entry['download_url'], headers={'User-Agent': 'PerturbCompass-reproduction'})
        with urllib.request.urlopen(request, timeout=120) as response, temporary.open('xb') as output:
            while block := response.read(16 * 1024 * 1024):
                output.write(block)
        if temporary.stat().st_size != entry['size_bytes'] or checksum(temporary) != entry['checksums']['sha256']:
            raise RuntimeError(f'Download failed byte checks; partial retained at {temporary}')
        temporary.rename(destination)
        print('Verified', entry['path'])


if __name__ == '__main__':
    main()
