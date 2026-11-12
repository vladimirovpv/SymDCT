#!/usr/bin/env python3
"""Check that sbom.cdx.json lists the same packages as a freshly exported SBOM.

Usage:
    uv export --format cyclonedx1.5 --output-file /tmp/sbom.new.json
    python scripts/check_sbom.py sbom.cdx.json /tmp/sbom.new.json

Only package names and versions are compared.  Time stamp, serial number and
tool version differ between two exports and are ignored.
"""

import json
import sys
from pathlib import Path


def components(path: str) -> set[tuple[str, str]]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return {(c["name"], c.get("version", "")) for c in data.get("components", [])}


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    committed, fresh = components(sys.argv[1]), components(sys.argv[2])
    if committed == fresh:
        print(f"SBOM is up to date ({len(committed)} components).")
        return 0
    for name, version in sorted(fresh - committed):
        print(f"missing in {sys.argv[1]}: {name} {version}")
    for name, version in sorted(committed - fresh):
        print(f"no longer in uv.lock: {name} {version}")
    print("Regenerate with: uv export --format cyclonedx1.5 --output-file sbom.cdx.json")
    return 1


if __name__ == "__main__":
    sys.exit(main())
