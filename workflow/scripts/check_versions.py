#!/usr/bin/env python3
"""Preflight: verify the active Python env satisfies the pins in requirements.txt.

The pipeline's Python steps (starsolo_to_h5ad, merge, downstream) import a small
subset of the project's packages. This reads the exact pins for that subset from
requirements.txt and compares them to what's installed, exiting non-zero on any
mismatch or missing package so the workflow stops before doing real work.
"""
import argparse
import importlib.metadata as md
import re
import sys

# Distribution names the pipeline actually depends on (import name may differ,
# e.g. pyyaml -> yaml). We compare by distribution name against requirements.txt.
REQUIRED = ["anndata", "scanpy", "numpy", "pandas", "scipy", "h5py",
            "pyyaml", "requests"]

PIN_RE = re.compile(r"^([A-Za-z0-9_.\-]+)==([0-9][^\s\\;]*)")


def parse_pins(requirements_path):
    pins = {}
    with open(requirements_path) as fh:
        for line in fh:
            m = PIN_RE.match(line.strip())
            if m:
                pins[m.group(1).lower()] = m.group(2)
    return pins


def installed_version(dist):
    for name in (dist, dist.replace("pyyaml", "PyYAML")):
        try:
            return md.version(name)
        except md.PackageNotFoundError:
            continue
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("requirements", help="path to requirements.txt")
    ap.add_argument("--warn-only", action="store_true",
                    help="report mismatches but exit 0")
    args = ap.parse_args()

    pins = parse_pins(args.requirements)
    problems = []
    print(f"{'package':12} {'required':12} {'installed':12} status")
    for pkg in REQUIRED:
        want = pins.get(pkg)
        have = installed_version(pkg)
        if want is None:
            status = "not pinned"          # not in requirements.txt; skip
        elif have is None:
            status = "MISSING"; problems.append(pkg)
        elif have != want:
            status = "MISMATCH"; problems.append(pkg)
        else:
            status = "ok"
        print(f"{pkg:12} {str(want):12} {str(have):12} {status}")

    if problems and not args.warn_only:
        print(f"\nEnvironment does not satisfy requirements.txt: {problems}",
              file=sys.stderr)
        sys.exit(1)
    print("\nEnvironment OK" if not problems else "\n(warn-only) mismatches ignored")


if __name__ == "__main__":
    main()
