# Reproducing PerturbCompass

Use the environment in `environment.yml` and an editable installation (`python -m pip install --no-deps -e .`). The original Python and numerical-library versions are recorded there because randomized SVD, clustering and linear algebra may vary across environments.

**Verification without raw inputs**

```bash
python -m pytest -q
python -m perturbcompass verify
python -m perturbcompass results
```

These commands validate retained evidence hashes and independently reconstruct primary regret point estimates and the 256-anchor identity from the reference CSVs. They do not refit the models or rebuild cell-level counts.

**Raw data**

Raw data are excluded from Git and release archives. `data/manifest.json` specifies each required file, its size and SHA-256. Existing local inputs have been renamed by study and their bytes retained. Download public inputs using `python scripts/download_data.py --list` and then `--dataset STUDY_NAME`, or provide the byte-identical author artifacts at the manifest paths. Download destinations are staged and must pass size and hash checks before acceptance.

```bash
python -m perturbcompass verify --check-local-data
python -m perturbcompass verify --include-raw
```

Full raw hashing reads tens of gigabytes and can take time. The package never silently substitutes a different artifact or regenerates the eligibility split from target outcomes.

**Model analyses**

```bash
python -m perturbcompass run cd4t_confirmation
python -m perturbcompass run replicate_consistency
python -m perturbcompass run source_ood
python -m perturbcompass run fallback_sensitivity
python -m perturbcompass run anchor_budget
```

Each study uses a descriptive configuration from `configs/`, fixed eligibility/partition files from `data/processed/`, and the same model hyperparameters, numerical seeds and inferential rules as the retained working-tree implementation. New outputs are written to `results/reproduced/STUDY_NAME/`. Reference results and historical seals are not overwritten.

The fixed source-cluster labels and farthest-first anchor order are provided with hashes. The generating pure functions are in `design.py`; the original generating source is preserved in provenance. This release does not re-execute historical freeze scripts or claim to recreate their dates.

Reorganization verification includes regression tests and a numerical comparison with the original five-family predictions on a deterministic synthetic fixture. See `docs/provenance/refactor_validation.json` for checks actually completed. Full raw-to-result execution is listed separately and must not be inferred from those checks.

**Creating a portable release**

```bash
python scripts/package_release.py --output /path/outside/repository/PerturbCompass_source.zip
```

The archive includes source, protocols, reference tables, partitions, documentation and provenance. It excludes `.git`, raw data, generated caches, reproduction outputs and private recovery locations. It does not push to GitHub. See `release_status.md` for author-controlled release metadata.
