# Provenance and reorganization

The repository was reorganized after the reported results were available. The organizational change does not constitute a new preregistration, a new external validation, or evidence that historical outcome access can be independently reconstructed from hashes alone.

`original_protocol.yaml`, amendments, numeric-access seals and reference results are retained byte-for-byte. They contain historical identifiers and old relative paths. Resolve those paths with `path_mapping.json`; do not edit the frozen documents to make their paths appear current.

`file_reorganization.csv` accounts for every original working-tree file. Exploratory documents and protocols are embedded in `../research_history.html`. Original code, large exploratory inputs and temporary artifacts remain recoverable in the private backup outside this repository.

`immutable_artifacts.json` binds retained original records, reference tables and fixed designs. `reorganized_code_manifest.json` binds the new package, tests and active configurations. These are distinct checks: the second manifest is explicitly post-result and does not replace the original freeze dates.

`source_snapshots/` preserves original versions of the paper-related computational source as non-executable text. The byte hashes connect to the historical seals. The active implementation lives in `src/perturbcompass/` with reusable functions separated from study entry points.

Run `python -m perturbcompass verify` to check the retained small artifacts and `python -m perturbcompass results` to reconcile the primary result tables. Reproduction checks selected raw inputs against their historical SHA-256 hashes before model fitting. See `../reproducibility.md`.
