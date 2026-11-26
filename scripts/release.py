#!/usr/bin/env python3
"""Prepare a release of SymDCT.

Usage
    uv run python scripts/release.py 0.1.2            prepare release 0.1.2
    uv run python scripts/release.py 0.1.2 --dry-run  show the plan only
    uv run python scripts/release.py --check [TAG]    check that versions agree

The release step does the following in this order.
    1. Checks the version string, the git tag and a clean working tree.
    2. Writes the version into pyproject.toml, CITATION.cff and README.md,
       and the release date into CITATION.cff.
    3. Runs "uv lock" (updates uv.lock) and exports sbom.cdx.json from uv.lock.
    4. Runs the test suite.
    5. Commits the changed files and creates the annotated tag vX.Y.Z.
Nothing is pushed. The Zenodo DOI of a version exists only after the release is
published on GitHub, so CITATION.cff and README.md carry the DOI of all versions.

The --check mode compares the versions in pyproject.toml, CITATION.cff, README.md
and uv.lock (and with TAG if given). It is meant for CI on a published release.
"""

from __future__ import annotations

import argparse
import datetime as dt
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PYPROJECT = ROOT / "pyproject.toml"
CITATION = ROOT / "CITATION.cff"
README = ROOT / "README.md"
LOCK = ROOT / "uv.lock"
SBOM = ROOT / "sbom.cdx.json"
RELEASE_FILES = [PYPROJECT, CITATION, README, LOCK, SBOM]

SEMVER = re.compile(r"^\d+\.\d+\.\d+$")
RE_PYPROJECT = re.compile(r'(?m)^version = "(\d+\.\d+\.\d+)"$')
RE_CITATION = re.compile(r"(?m)^version: (\d+\.\d+\.\d+)$")
RE_DATE = re.compile(r"(?m)^date-released: (\S+)$")
RE_README = re.compile(r"version (\d+\.\d+\.\d+),")
RE_LOCK = re.compile(r'name = "symdct"\nversion = "(\d+\.\d+\.\d+)"')


def run(cmd: list[str], check: bool = True) -> subprocess.CompletedProcess:
    print("$", " ".join(cmd))
    return subprocess.run(cmd, cwd=ROOT, check=check, text=True)


def capture(cmd: list[str]) -> str:
    return subprocess.run(cmd, cwd=ROOT, check=True, text=True, capture_output=True).stdout


def find(regex: re.Pattern, path: Path) -> str:
    match = regex.search(path.read_text(encoding="utf-8"))
    return match.group(1) if match else "<not found>"


def versions() -> dict[str, str]:
    return {
        "pyproject.toml": find(RE_PYPROJECT, PYPROJECT),
        "CITATION.cff": find(RE_CITATION, CITATION),
        "README.md": find(RE_README, README),
        "uv.lock": find(RE_LOCK, LOCK),
    }


def as_tuple(version: str) -> tuple[int, ...]:
    return tuple(int(part) for part in version.split("."))


def check(tag: str | None) -> int:
    found = versions()
    if tag:
        found["tag"] = tag.removeprefix("v")
    for name, value in found.items():
        print(f"{name:16s} {value}")
    if len(set(found.values())) == 1 and "<not found>" not in found.values():
        print("Versions agree.")
        return 0
    print("Versions differ.")
    return 1


def set_group(path: Path, regex: re.Pattern, value: str) -> None:
    """Replace group 1 of the single match of ``regex`` in ``path``."""
    text = path.read_text(encoding="utf-8")
    matches = list(regex.finditer(text))
    if len(matches) != 1:
        raise SystemExit(f"Expected exactly one match for {regex.pattern!r} in {path.name}")
    m = matches[0]
    path.write_text(text[: m.start(1)] + value + text[m.end(1) :], encoding="utf-8")


def release(version: str, date: str, dry_run: bool, skip_tests: bool, allow_dirty: bool) -> int:
    if not SEMVER.match(version):
        raise SystemExit(f"Version must have the form X.Y.Z, got {version!r}")
    current = find(RE_PYPROJECT, PYPROJECT)
    if as_tuple(version) <= as_tuple(current):
        raise SystemExit(f"New version {version} must be larger than the current version {current}")
    tag = f"v{version}"
    if capture(["git", "tag", "--list", tag]).strip():
        raise SystemExit(f"Tag {tag} exists")
    dirty = capture(["git", "status", "--porcelain", "--untracked-files=no", "--", "."]).strip()
    if dirty and not allow_dirty:
        raise SystemExit("Tracked files have uncommitted changes. Commit them first.\n" + dirty)

    print(f"Release {current} -> {version}, date {date}")
    print("Files:", ", ".join(p.name for p in RELEASE_FILES))
    if dry_run:
        print("Dry run, nothing written.")
        return 0

    set_group(PYPROJECT, RE_PYPROJECT, version)
    set_group(CITATION, RE_CITATION, version)
    set_group(CITATION, RE_DATE, date)
    set_group(README, RE_README, version)

    run(["uv", "lock"])
    run(["uv", "export", "--locked", "--format", "cyclonedx1.5", "--output-file", str(SBOM)])
    if check(tag) != 0:
        raise SystemExit("Versions disagree after the update.")
    if not skip_tests:
        run(["uv", "run", "pytest", "-q"])

    run(["git", "add", *[str(p) for p in RELEASE_FILES]])
    run(["git", "commit", "-m", f"Release {tag}"])
    run(["git", "tag", "-a", tag, "-m", f"SymDCT {version}"])
    print(f"\nPrepared {tag}. Next steps")
    print("  1. git push origin <branch> && git push origin", tag)
    print("  2. Mirror public/ to the SymDCT repository and publish the GitHub release.")
    print("  3. Note the Zenodo DOI of the version for papers and the release notes.")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("version", nargs="?", help="new version X.Y.Z")
    parser.add_argument("--check", nargs="?", const="", metavar="TAG", help="check that versions agree")
    parser.add_argument("--date", default=dt.date.today().isoformat(), help="release date, default today")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--skip-tests", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()
    if args.check is not None:
        return check(args.check or None)
    if not args.version:
        parser.error("version is required")
    return release(args.version, args.date, args.dry_run, args.skip_tests, args.allow_dirty)


if __name__ == "__main__":
    sys.exit(main())
