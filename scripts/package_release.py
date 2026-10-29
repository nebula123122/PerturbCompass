"""Package paper code and evidence while preserving portable archive paths."""
from __future__ import annotations

import argparse
from pathlib import Path, PurePath
import zipfile

ROOT = Path(__file__).resolve().parents[1]
DIRECTORIES = {'.github', 'src', 'tests', 'configs', 'docs', 'manuscript', 'figures', 'scripts'}
ROOT_FILES = {
    'README.md', 'README.zh-CN.md', '.gitignore', '.gitattributes', 'pyproject.toml',
    'environment.yml', 'LICENSE', 'LICENSE.md', 'LICENSE.txt', 'COPYING',
    'CITATION.cff', 'CITATION.md', 'CONTRIBUTING.md',
}


def is_release_path(relative: PurePath) -> bool:
    if not relative.parts or relative.is_absolute() or '..' in relative.parts:
        return False
    if any(part in {'__pycache__', '.pytest_cache', '.git', '.venv'} or part.endswith('.egg-info') for part in relative.parts):
        return False
    name = relative.as_posix()
    return (relative.parts[0] in DIRECTORIES or name in ROOT_FILES or name == 'data/manifest.json'
            or name.startswith(('data/processed/', 'results/reference/')))


def release_files():
    for path in sorted(ROOT.rglob('*')):
        if not path.is_file():
            continue
        relative = path.relative_to(ROOT)
        if is_release_path(relative):
            if path.is_symlink():
                raise RuntimeError(f'Release cannot contain symbolic links: {relative}')
            yield path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.is_relative_to(ROOT):
        parser.error('Choose an output path outside the repository.')
    if output.exists():
        parser.error('Output already exists; choose a new filename.')
    output.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with zipfile.ZipFile(output, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for path in release_files():
            archive.write(path, 'PerturbCompass/' + path.relative_to(ROOT).as_posix())
            count += 1
    print(f'{output}: {count} files; {output.stat().st_size} bytes')


if __name__ == '__main__':
    main()
