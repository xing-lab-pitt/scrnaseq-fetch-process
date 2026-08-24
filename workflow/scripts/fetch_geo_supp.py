#!/usr/bin/env python3
"""List or download the SUPPLEMENTARY files a GEO submitter attached to a series
or sample — the author-processed `.h5ad`, `.rds`, `barcodes/features/matrix`,
`.loom` and so on that live next to the record on the GEO FTP mirror.

This is deliberately NOT part of the Snakemake DAG. The pipeline's own product
is a `.h5ad` it builds itself from raw reads (FASTQ -> STARsolo -> h5ad), with a
known reference, a known chemistry and velocity layers. A supplementary h5ad is
whatever the authors happened to produce: different reference, different
filtering, usually no spliced/unspliced layers. The two are not interchangeable,
so this script keeps the author-made files on a separate, explicit path.

    # what is attached? (default -- lists, downloads nothing)
    python fetch_geo_supp.py GSE218661

    # fetch just the h5ad files
    python fetch_geo_supp.py GSE218661 -o results/geo_supp --download \
        --pattern '*.h5ad.gz' --gunzip

Downloads resume (`curl -C -`) and are size-verified against the server's
Content-Length, so an interrupted multi-GB transfer is restartable and a
truncated file is never silently accepted. A file already present at the right
size is skipped.

Exit status: 0 all requested files present and verified, 1 a transfer failed,
2 the accession or pattern matched nothing.
"""

import argparse
import fnmatch
import gzip
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo"
TIMEOUT = 60

# Which FTP subtree an accession lives under.
_SUBTREE = {"GSE": "series", "GSM": "samples", "GPL": "platforms"}


def suppl_url(accession: str) -> str:
    """Map an accession to its supplementary directory URL.

    GEO buckets records by masking the last three digits, so GSE218661 lands in
    `series/GSE218nnn/GSE218661/suppl/`.
    """
    acc = accession.strip().upper()
    m = re.fullmatch(r"(GSE|GSM|GPL)(\d+)", acc)
    if not m:
        raise ValueError(f"not a GEO accession: {accession!r} (expected GSE/GSM/GPL + digits)")
    prefix, digits = m.groups()
    # Accessions below 1000 keep an empty stem, so GSE1 lives in `GSEnnn/GSE1/`.
    bucket = f"{prefix}{digits[:-3]}nnn"
    return f"{FTP_BASE}/{_SUBTREE[prefix]}/{bucket}/{acc}/suppl/"


def list_files(dir_url: str) -> list:
    """Return the file names in an Apache directory index, parents excluded."""
    try:
        with urllib.request.urlopen(dir_url, timeout=TIMEOUT) as r:
            html = r.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"{dir_url} -> HTTP {e.code} (no supplementary files, or wrong accession)")
    names = []
    for href in re.findall(r'href="([^"?/][^"]*)"', html):
        if href.startswith(("http:", "https:", "#")):
            continue
        name = urllib.parse.unquote(href)
        if name not in names:
            names.append(name)
    return names


def remote_size(url: str):
    """Exact byte size from a HEAD request, or None if the server won't say."""
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
            n = r.headers.get("Content-Length")
            return int(n) if n else None
    except urllib.error.HTTPError:
        return None


def human(n) -> str:
    if n is None:
        return "?"
    for unit in ("B", "K", "M", "G", "T"):
        if n < 1024 or unit == "T":
            return f"{n:.0f}{unit}" if unit == "B" else f"{n:.1f}{unit}"
        n /= 1024.0


def fetch_one(url: str, dest: Path, size, quiet: bool = False) -> bool:
    """Resume-capable download plus a size check. True when dest is complete."""
    if size is not None and dest.exists() and dest.stat().st_size == size:
        print(f"[skip] {dest.name} already complete ({human(size)})")
        return True
    dest.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["curl", "-fL", "-C", "-", "-o", str(dest), url]
    if quiet:
        cmd.insert(1, "-sS")
    print(f"[get ] {dest.name} ({human(size)})")
    rc = subprocess.call(cmd)
    if rc != 0:
        print(f"[FAIL] curl exit {rc} on {url}", file=sys.stderr)
        return False
    if size is not None and dest.stat().st_size != size:
        print(f"[FAIL] {dest.name}: got {dest.stat().st_size} bytes, expected {size}", file=sys.stderr)
        return False
    return True


def gunzip(path: Path) -> Path:
    """Decompress in place, keeping the .gz until the plain file is complete."""
    out = path.with_suffix("")
    if out.exists():
        print(f"[skip] {out.name} already decompressed")
        return out
    print(f"[gzip] {path.name} -> {out.name}")
    tmp = out.with_suffix(out.suffix + ".part")
    with gzip.open(path, "rb") as fi, open(tmp, "wb") as fo:
        shutil.copyfileobj(fi, fo, length=1 << 22)
    tmp.rename(out)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("accession", help="GSE / GSM / GPL accession")
    ap.add_argument("-o", "--outdir", default=".", help="download destination (default: cwd)")
    ap.add_argument("--download", action="store_true",
                    help="actually transfer; without it the script only lists")
    ap.add_argument("--pattern", default="*",
                    help="glob over file names, e.g. '*.h5ad.gz' (default: all)")
    ap.add_argument("--gunzip", action="store_true", help="decompress .gz after a verified download")
    ap.add_argument("--quiet", action="store_true", help="suppress curl's progress meter")
    args = ap.parse_args()

    try:
        base = suppl_url(args.accession)
        names = list_files(base)
    except (ValueError, RuntimeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    picked = [n for n in names if fnmatch.fnmatch(n, args.pattern)]
    if not picked:
        print(f"error: {args.accession}: {len(names)} supplementary file(s), none matching "
              f"{args.pattern!r}", file=sys.stderr)
        for n in names:
            print(f"  {n}", file=sys.stderr)
        return 2

    sizes = {n: remote_size(base + n) for n in picked}
    print(f"{args.accession} supplementary files at {base}")
    for n in picked:
        print(f"  {human(sizes[n]):>8}  {n}")
    total = sum(s for s in sizes.values() if s)
    print(f"  {'-' * 8}")
    print(f"  {human(total):>8}  total ({len(picked)} file(s))")

    if not args.download:
        print("\n(listing only -- pass --download to transfer)")
        return 0

    outdir = Path(args.outdir).expanduser()
    failed = []
    for n in picked:
        dest = outdir / n
        if fetch_one(base + n, dest, sizes[n], quiet=args.quiet):
            if args.gunzip and dest.suffix == ".gz":
                gunzip(dest)
        else:
            failed.append(n)

    if failed:
        print(f"\n{len(failed)} of {len(picked)} failed: {', '.join(failed)}", file=sys.stderr)
        return 1
    print(f"\nOK: {len(picked)} file(s) in {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
