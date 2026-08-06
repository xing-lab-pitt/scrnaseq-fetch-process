#!/usr/bin/env python3
"""Resolve the STAR index the pipeline will align against — as a pipeline OUTPUT.

Principle: the index and the aligner must come from the same environment. STAR
refuses to load an index whose on-disk `versionGenome` differs from what the
running STAR binary speaks (e.g. STAR 2.7.10a reads a 2.7.4a index fine, but a
2.7.1a index aborts with "Genome version ... is INCOMPATIBLE"). The compatibility
key is the index's `versionGenome`, NOT the binary's version string — two STAR
releases can share one genome version.

So we never let `starsolo` consume a hand-supplied index directly. Instead this
step produces `WORK/star_index`, which the pipeline vouches for:

  * supplied index whose versionGenome MATCHES this STAR  -> symlink it in (instant)
  * supplied index that MISMATCHES, or no index supplied  -> build from fasta+gtf

We discover "what this STAR speaks" empirically by building a throwaway probe
index (~0.03s) rather than hardcoding a version-compatibility table.
"""
import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def _star() -> str:
    exe = shutil.which("STAR")
    if not exe:
        sys.exit("STAR not found on PATH — cannot resolve the index.")
    return exe


def _symlink(out: Path, target: Path):
    """Point the pipeline-owned `out` at an index dir that lives elsewhere."""
    if out.is_symlink() or out.exists():
        if out.is_symlink() or out.is_file():
            out.unlink()
        else:
            shutil.rmtree(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.symlink_to(target.resolve())


def _genome_version(index_dir: Path):
    """Read `versionGenome` from an index's genomeParameters.txt, or None."""
    params = index_dir / "genomeParameters.txt"
    if not params.exists():
        return None
    for line in params.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "versionGenome":
            return parts[1]
    return None


def probe_genome_version() -> str:
    """Build a tiny throwaway index to learn this STAR's native genomeVersion."""
    with tempfile.TemporaryDirectory(prefix="_starprobe.") as tmp:
        tmp = Path(tmp)
        (tmp / "tiny.fa").write_text(">probe\n" + "ACGT" * 500 + "\n")
        idx = tmp / "idx"
        idx.mkdir()
        subprocess.run(
            [_star(), "--runMode", "genomeGenerate", "--genomeDir", str(idx),
             "--genomeFastaFiles", str(tmp / "tiny.fa"),
             "--genomeSAindexNbases", "4", "--runThreadN", "2"],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, cwd=tmp)
        ver = _genome_version(idx)
        if ver is None:
            sys.exit("probe index built but had no versionGenome — cannot resolve.")
        return ver


def build_index(out: Path, fasta: str, gtf: str, overhang, threads):
    for label, p in (("fasta", fasta), ("gtf", gtf)):
        if not p or not Path(p).exists():
            sys.exit(f"Need to build the STAR index but reference.{label} is "
                     f"missing ({p!r}); set it in config or supply a compatible index.")
    out.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [_star(), "--runMode", "genomeGenerate", "--runThreadN", str(threads),
         "--genomeDir", str(out), "--genomeFastaFiles", fasta,
         "--sjdbGTFfile", gtf, "--sjdbOverhang", str(overhang)],
        check=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--out", required=True, help="index dir the pipeline will use")
    ap.add_argument("--supplied", default="", help="optional prebuilt index path from config")
    ap.add_argument("--fasta", default="")
    ap.add_argument("--gtf", default="")
    ap.add_argument("--overhang", default="100")
    ap.add_argument("--threads", default="8")
    args = ap.parse_args()

    out = Path(args.out)
    native = probe_genome_version()
    print(f"[resolve_star_index] this STAR speaks genomeVersion {native}")

    supplied = Path(args.supplied) if args.supplied else None
    if supplied and supplied != out and _genome_version(supplied) is not None:
        have = _genome_version(supplied)
        if have == native:
            # Compatible -> reuse in place via symlink (no copy, no rebuild).
            _symlink(out, supplied)
            print(f"[resolve_star_index] supplied index {supplied} is genomeVersion "
                  f"{have} == native; symlinked -> {out} (no rebuild).")
            return
        print(f"[resolve_star_index] supplied index {supplied} is genomeVersion "
              f"{have}, INCOMPATIBLE with this STAR ({native}); rebuilding from "
              f"fasta+gtf so the index and aligner match.")
    elif supplied and supplied != out:
        print(f"[resolve_star_index] supplied index {supplied} not found / not a "
              f"STAR index; building from fasta+gtf.")
    else:
        print("[resolve_star_index] no prebuilt index supplied; building from fasta+gtf.")

    # Where the index bytes physically land. If a --supplied path was given, build
    # there so the index is DURABLE (survives the ephemeral workdir and is reused on
    # the next run — the compatible branch above will then just symlink it, no
    # rebuild). WORK/star_index is then a symlink to it, so the pipeline still owns
    # the path the aligner reads. With no --supplied path, build directly into out.
    dest = supplied if (supplied and supplied != out) else out
    build_index(dest, args.fasta, args.gtf, args.overhang, args.threads)
    built = _genome_version(dest)
    if dest != out:
        _symlink(out, dest)
        print(f"[resolve_star_index] built index at {dest} (genomeVersion {built}); "
              f"symlinked -> {out}.")
    else:
        print(f"[resolve_star_index] built index at {out} (genomeVersion {built}).")


if __name__ == "__main__":
    main()
