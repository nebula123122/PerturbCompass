"""The release must contain evidence and notices on Windows and Unix."""
import importlib.util
from pathlib import Path, PurePosixPath, PureWindowsPath
import sys
import zipfile

from perturbcompass.runtime import ROOT

spec = importlib.util.spec_from_file_location('paper_release', ROOT / 'scripts/package_release.py')
release = importlib.util.module_from_spec(spec)
spec.loader.exec_module(release)


def test_portable_paths_include_designs_and_reference_tables():
    for cls in [PurePosixPath, PureWindowsPath]:
        for name in ['data/manifest.json', 'data/processed/study/partition.csv', 'results/reference/study/metrics.csv', '.github/workflows/tests.yml', '.gitattributes', 'LICENSE', 'CITATION.cff']:
            assert release.is_release_path(cls(name)), name


def test_private_data_and_generated_files_are_excluded():
    for cls in [PurePosixPath, PureWindowsPath]:
        for name in ['data/raw/counts.h5ad', '.git/config', 'results/reproduced/run.csv', 'src/__pycache__/module.pyc', 'src/example.egg-info/PKG-INFO', '../outside.txt']:
            assert not release.is_release_path(cls(name)), name


def test_archive_keeps_nested_paths_and_future_license(tmp_path, monkeypatch):
    project = tmp_path / 'project'; project.mkdir()
    for name in ['README.md', '.gitattributes', 'LICENSE', 'CITATION.cff', 'data/manifest.json', 'data/processed/study/partition.csv', 'results/reference/study/metrics.csv', 'data/raw/excluded.bin']:
        file = project / name; file.parent.mkdir(parents=True, exist_ok=True); file.write_text(name, encoding='utf-8')
    monkeypatch.setattr(release, 'ROOT', project)
    output = tmp_path / 'release.zip'
    monkeypatch.setattr(sys, 'argv', ['package_release.py', '--output', str(output)])
    release.main()
    with zipfile.ZipFile(output) as bundle:
        names = set(bundle.namelist())
        assert 'PerturbCompass/LICENSE' in names
        assert 'PerturbCompass/data/manifest.json' in names
        assert 'PerturbCompass/results/reference/study/metrics.csv' in names
        assert not any('\\' in name or '/raw/' in name for name in names)
