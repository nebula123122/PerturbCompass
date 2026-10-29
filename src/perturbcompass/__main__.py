"""Command-line entry points for paper verification and reproduction."""
from __future__ import annotations

import argparse
import importlib
import json
import sys
from pathlib import Path

from perturbcompass.runtime import ROOT, verify_raw_data, verify_small_artifacts

ANALYSES = [
    'k562_cross_study', 'k562_screen_validation', 'feng_external_validation',
    'feng_robustness', 'nadig_selector_validation', 'cd4t_confirmation',
    'replicate_consistency', 'source_ood', 'fallback_sensitivity', 'anchor_budget',
]


def main():
    parser = argparse.ArgumentParser(prog='perturbcompass')
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('analyses', help='List paper analyses')
    verify = sub.add_parser('verify', help='Verify code, frozen records, partitions and reference tables')
    verify.add_argument('--include-raw', action='store_true', help='Also hash all local raw inputs')
    verify.add_argument('--check-local-data', action='store_true', help='Check local raw input sizes')
    sub.add_parser('results', help='Recompute primary conclusions from immutable reference tables')
    run = sub.add_parser('run', help='Reproduce an analysis; writes only results/reproduced/')
    run.add_argument('analysis', choices=ANALYSES)
    args = parser.parse_args()
    if args.command == 'analyses':
        print('\n'.join(ANALYSES))
    elif args.command == 'verify':
        result = {'small_artifacts': verify_small_artifacts()}
        if args.include_raw or args.check_local_data:
            result['raw_data'] = verify_raw_data(include_hashes=args.include_raw)
        print(json.dumps(result, indent=2))
    elif args.command == 'results':
        from perturbcompass.reference_results import audit_reference_results
        print(json.dumps(audit_reference_results(), indent=2))
    else:
        verify_small_artifacts()
        output = ROOT / 'results/reproduced' / args.analysis
        output.mkdir(parents=True, exist_ok=True)
        # Retained analyses inherited argparse. Supply only their default config.
        sys.argv = [args.analysis]
        importlib.import_module('perturbcompass.analyses.' + args.analysis).main()


if __name__ == '__main__':
    main()
