#!/usr/bin/env python3
"""
Deterministically partition a BAM into N pseudo-replicates by hashing read QNAME.
Paired-end safe: R1 and R2 share a QNAME, so both mates land in the same bucket.

Usage:
    split_bam_pseudoreps.py --input INPUT.bam --output-prefix PFX --n-buckets 3 --seed 42

Outputs:
    PFX_pseudo1.bam, PFX_pseudo2.bam, PFX_pseudo3.bam (and indexes)
"""
import argparse
import hashlib
import os
import sys
from typing import List


def qname_to_bucket(qname: str, n_buckets: int, seed: int) -> int:
    """Map a QNAME to a bucket index in [0, n_buckets).

    Deterministic for a given (qname, n_buckets, seed). Uses MD5 of the
    seed-prefixed qname; the first 8 bytes are reinterpreted as a 64-bit int
    and reduced modulo n_buckets.
    """
    payload = f"{seed}:{qname}".encode("utf-8")
    digest = hashlib.md5(payload).digest()
    val = int.from_bytes(digest[:8], byteorder="big", signed=False)
    return val % n_buckets


def split_bam(input_bam: str, output_prefix: str, n_buckets: int, seed: int) -> List[str]:
    """Partition BAM records into n_buckets output BAMs by qname hash. Returns output paths."""
    import pysam

    in_bam = pysam.AlignmentFile(input_bam, "rb")
    out_paths = [f"{output_prefix}_pseudo{i+1}.bam" for i in range(n_buckets)]
    out_bams = [pysam.AlignmentFile(p, "wb", template=in_bam) for p in out_paths]

    try:
        for rec in in_bam.fetch(until_eof=True):
            bucket = qname_to_bucket(rec.query_name, n_buckets, seed)
            out_bams[bucket].write(rec)
    finally:
        in_bam.close()
        for ob in out_bams:
            ob.close()

    # Index each output
    for p in out_paths:
        pysam.index(p)

    return out_paths


def main():
    parser = argparse.ArgumentParser(description="Deterministic BAM pseudo-replicate splitter")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--n-buckets", type=int, default=3)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    if args.n_buckets < 2:
        print("ERROR: --n-buckets must be at least 2", file=sys.stderr)
        sys.exit(1)

    paths = split_bam(args.input, args.output_prefix, args.n_buckets, args.seed)
    print(f"Wrote {len(paths)} pseudo-replicate BAMs:")
    for p in paths:
        print(f"  {p}")


if __name__ == "__main__":
    main()
