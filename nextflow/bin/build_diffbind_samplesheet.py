#!/usr/bin/env python3
"""Fan a cohort's differential-peaks rows out into one DiffBind samplesheet per case."""
import argparse
import csv
import json
import re
import sys

_PSEUDO_SUFFIX = re.compile(r"_pseudo\d+$")
RESERVED_LABELS = {"control", "background"}


def group_rows(rows):
    """Split diff rows into (control_rows, {case_label: [rows]})"""
    controls = [r for r in rows if r[1] == "control"]
    cases = {}
    for r in rows:
        if r[1] in RESERVED_LABELS:
            continue
        cases.setdefault(r[1], []).append(r)
    return controls, cases

def format_case_summary(cases):
    """Print to case plan to stdout"""
    labels = ["control"] + sorted(cases)
    lines = [f"{len(labels)} distinct treatment labels found: {', '.join(labels)}",
             "", "Cases:"]
    for i, case in enumerate(sorted(cases), start=1):
        rows = cases[case]
        pseudo = sum(1 for r in rows if _PSEUDO_SUFFIX.search(r[0]))
        real = len(rows) - pseudo
        lines.append(f"  [{i}] control vs {case} treatment ({real} real, {pseudo} pseudo reps)")
    return "\n".join(lines)


def build_case_samplesheet(controls, case_rows, peaks_map):
    """Pipeline-samplesheet rows for one case: shared control + treatment"""
    header = ["sample_id", "condition", "replicate", "bam_path", "peaks_path"]
    out = [header]
    for cond_label, rows in (("control", controls), ("treatment", case_rows)):
        for sid, _orig_cond, rep, bam in rows:
            parent = _PSEUDO_SUFFIX.sub("", sid)   # pseudo-reps reuse the parent's peaks
            peaks = peaks_map.get(sid) or peaks_map.get(parent) or ""
            out.append([sid, cond_label, rep, bam, peaks])
    return out

def _write_csv(path, table):
    with open(path, "w", newline="") as fh:
        csv.writer(fh).writerows(table)

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", required=True, help="JSON list of [id, condition, replicate, bam]")
    ap.add_argument("--peaks", required=True, help="JSON list of [id, peaks_path] pairs")
    ap.add_argument("--output-dir", default=".")
    args = ap.parse_args(argv)

    with open(args.rows) as fh:
        rows = json.load(fh)
    with open(args.peaks) as fh:
        peaks_map = {p[0]: p[1] for p in json.load(fh)}

    controls, cases = group_rows(rows)

    if not controls or not cases:
        reason = "no control row found" if not controls else "no treatment cases found"
        labels = sorted({r[1] for r in rows})
        sys.stderr.write(f"ERROR: invalid differential-peaks design: {reason}\n")
        sys.stderr.write(f"Diff-stage condition labels seen: {', '.join(labels)}\n")
        return 1

    print(format_case_summary(cases))

    for case in sorted(cases):
        table = build_case_samplesheet(controls, cases[case], peaks_map)
        _write_csv(f"{args.output_dir}/diffbind_samplesheet_{case}.csv", table)
    return 0


if __name__ == "__main__":
    sys.exit(main())
