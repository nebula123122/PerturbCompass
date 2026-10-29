"""Create a portable source release without raw data or private Git history."""
from __future__ import annotations

import argparse
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
DIRECTORIES = {'src', 'tests', 'configs', 'docs', 'manuscript', 'figures', 'scripts'}
ROOT_FILES = {'README.md', 'README.zh-CN.md', '.gitignore', 'pyproject.toml', 'environment.yml'}


def release_files():
    for path in sorted(ROOT.rglob('*')):
        if not path.is_file():
            continue
        relative = path.relative_to(ROOT)
        if any(part in {'__pycache__', '.pytest_cache', '.git', '.venv'} or part.endswith('.egg-info') for part in relative.parts):
            continue
        if relative.parts[0] in DIRECTORIES or str(relative) in ROOT_FILES or str(relative) == 'data/manifest.json' or str(relative).startswith(('data/processed/', 'results/reference/')):
            if path.is_symlink():
                raise RuntimeError(f'Release cannot contain symbolic links: {relative}')
            yield path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.is_relative_to(ROOT):
        parser.error('Choose an output path outside the repository to avoid packaging an archive recursively.')
    if output.exists():
        parser.error('Output already exists; choose a new filename.')
    output.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with zipfile.ZipFile(output, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for path in release_files():
            archive.write(path, 'PerturbCompass/' + str(path.relative_to(ROOT)))
            count += 1
    print(f'{output}: {count} files; {output.stat().st_size} bytes')


if __name__ == '__main__':
    main()
