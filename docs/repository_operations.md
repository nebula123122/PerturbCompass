# Server and Git workflow

The canonical research checkout is `PerturbCompass/` (the primary server working directory). It is on `main` and tracks `origin/main` at [nebula123122/PerturbCompass](https://github.com/nebula123122/PerturbCompass). Subsequent edits, commits and pushes should use this checkout.

The initial upload used a separate clean checkout to avoid publishing exploratory Git history. The canonical checkout was subsequently connected to that published history after all tracked file contents were compared. The former local Git metadata and the initial upload checkout remain recoverable outside the project. Their private storage paths and authentication keys are excluded from the release.

**Inspect and update**

```bash
cd /path/to/PerturbCompass
git status
git pull --ff-only
```

After making reviewed changes, refresh the current code manifest when code, configuration, tests, environment or CI changes. Then run the checks:

```bash
python scripts/refresh_integrity_manifest.py --reason "Describe the reviewed change"
python -m pytest -q
python -m perturbcompass verify
python -m perturbcompass results
git diff
git add --all
git commit -m "Describe the completed change"
git push origin main
```

Raw inputs, reproduction outputs and private recovery records stay outside Git tracking. Frozen historical records are preserved byte-for-byte. A new post-result method experiment requires its own scientific classification; a maintenance commit or refreshed hash is not a new preregistration.

Version tags mark author-approved milestones. The package metadata currently says `0.1.0`; no release tag is created automatically by these commands. The code license remains deferred by the author.
