"""
Synchronizer CLI

Usage:
  - as a module from the project "public" folder:
	  python -m src.common_utils.synchronizer /path/to/energetics_template.yaml /path/to/production_energetics.yaml

  - or directly (when PYTHONPATH includes the project `public` directory):
	  python src/common_utils/synchronizer.py /path/to/energetics_template.yaml /path/to/production_energetics.yaml

This script calls `synchronize_energetics(template, production)` from
`src.common_utils.energetics_sync` and writes a new production YAML that is
aligned with the template. It never overwrites the original production file.

Options:
  --out-dir DIR        Directory where the synchronized file will be written.
					   If omitted, the new file is written next to the
					   original production file.
  --length-tolerance Tolerance (float) passed to the synchroniser to match
					   jump lengths. If omitted, the default from the helper
					   is used.
  --dry-run            Do not write a file; show what would be written.
  -v / --verbose       Print more information.

Example:
  python -m src.common_utils.synchronizer examples/hcp_vac/output/energetics_template.yaml \
	  examples/hcp_vac/input/energetics.yaml --out-dir examples/hcp_vac/input
"""

from __future__ import annotations

import argparse
import json
import pprint
import sys
from typing import Optional

from src.common_utils.energetics_io import synchronize_energetics


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
	p = argparse.ArgumentParser(
		description="Synchronize a production energetics YAML with a template (do not overwrite originals)."
	)
	p.add_argument("template", help="Path to energetics_template.yaml")
	p.add_argument("production", help="Path to production energetics YAML to synchronize")
	p.add_argument("--out-dir", help="Directory to write the synchronized file to (defaults to production file dir)")
	p.add_argument("--length-tolerance", type=float, help="Length tolerance to use when matching jumps")
	p.add_argument("--dry-run", action="store_true", help="Do not write the output file; only show the planned actions")
	p.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
	return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
	args = _parse_args(argv)

	if args.verbose:
		print("Template:", args.template)
		print("Production:", args.production)
		if args.out_dir:
			print("Out dir:", args.out_dir)
		if args.length_tolerance is not None:
			print("Length tolerance:", args.length_tolerance)

	try:
		new_path, report = synchronize_energetics(
			args.template, args.production, out_dir=args.out_dir, length_tolerance=args.length_tolerance
		)
	except Exception as exc:  # pragma: no cover - top-level CLI error handling
		print("Error while synchronizing energetics:\n", exc, file=sys.stderr)
		return 2

	if args.dry_run:
		print("Dry-run: no file written.")
		print("Planned output path:", new_path)
	else:
		print("Wrote:", new_path)

	# Nicely pretty-print the report (dict)
	print("Report:")
	pprint.pprint(report)

	# Also emit machine-readable JSON to stdout if verbose requested
	if args.verbose:
		try:
			print(json.dumps({"new_file": new_path, "report": report}, indent=2))
		except TypeError:
			# Some values in report may not be JSON serializable; ignore
			pass

	return 0


if __name__ == "__main__":
	raise SystemExit(main())
