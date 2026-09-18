#!/usr/bin/env python3
"""
Validator for DSSAT BIOCHAR.OUT output files.

Usage:
  python check_biochar_out.py BIOCHAR.OUT          # validate a real output
  python check_biochar_out.py --self-test          # run built-in unit tests
"""

import sys
import os
import io

# Expected column header tokens (subset — order and count must match)
EXPECTED_HEADER_TOKENS = [
    "@YEAR", "DOY", "DAS", "BIOCHC", "BCLBL", "BCSTB",
    "BIOCHN", "BCMINC", "DPH_L1", "SORBNH4", "NAPBIO", "CUMBCC", "CUMBCN",
]

# Column index in the data rows (0-based after splitting)
COL = {
    "YEAR":    0,
    "DOY":     1,
    "DAS":     2,
    "BIOCHC":  3,
    "BCLBL":   4,
    "BCSTB":   5,
    "BIOCHN":  6,
    "BCMINC":  7,
    "DPH_L1":  8,
    "SORBNH4": 9,
    "NAPBIO":  10,
    "CUMBCC":  11,
    "CUMBCN":  12,
}

ALLOWED_RANGE = {
    "BIOCHC":  (0.0, 1e7),
    "BCLBL":   (0.0, 1e7),
    "BCSTB":   (0.0, 1e7),
    "BIOCHN":  (0.0, 1e6),
    "BCMINC":  (0.0, 1e5),
    "DPH_L1":  (0.0, 2.0),
    "SORBNH4": (0.0, 1e5),
    "NAPBIO":  (0.0, 9999),
    "CUMBCC":  (0.0, 1e8),
    "CUMBCN":  (0.0, 1e7),
}


def check_file(path):
    errors = []
    warnings = []

    if not os.path.isfile(path):
        return [f"File not found: {path}"], []

    with open(path) as fh:
        lines = fh.readlines()

    if not lines:
        return ["File is empty"], []

    # First line: title
    if not lines[0].strip().startswith("*BIOCHAR"):
        warnings.append(f"Line 1 expected '*BIOCHAR...' title, got: {lines[0].rstrip()!r}")

    # Find header line (@YEAR ...)
    header_line = None
    header_idx = None
    for i, ln in enumerate(lines):
        if ln.strip().startswith("@YEAR"):
            header_line = ln
            header_idx = i
            break

    if header_line is None:
        return errors + ["No header line (@YEAR ...) found in file"], warnings

    tokens = header_line.split()
    missing = [t for t in EXPECTED_HEADER_TOKENS if t not in tokens]
    if missing:
        errors.append(f"Header missing columns: {missing}")

    # Validate data rows
    data_rows = []
    prev_biochc = None
    for i, ln in enumerate(lines[header_idx + 1:], start=header_idx + 2):
        ln = ln.rstrip()
        if not ln or ln.startswith("!") or ln.startswith("*") or ln.startswith("@"):
            continue
        parts = ln.split()
        if len(parts) < len(COL):
            warnings.append(f"Line {i}: too few columns ({len(parts)} < {len(COL)})")
            continue

        row = {}
        try:
            for name, idx in COL.items():
                row[name] = float(parts[idx])
        except (ValueError, IndexError) as exc:
            errors.append(f"Line {i}: parse error: {exc}")
            continue

        for col, (lo, hi) in ALLOWED_RANGE.items():
            val = row[col]
            if not (lo <= val <= hi):
                errors.append(
                    f"Line {i}: {col}={val:.4g} outside expected range [{lo}, {hi}]"
                )

        # BCLBL + BCSTB == BIOCHC (within tolerance)
        total = row["BCLBL"] + row["BCSTB"]
        if abs(total - row["BIOCHC"]) > max(1e-3 * max(total, 1.0), 0.01):
            errors.append(
                f"Line {i}: BCLBL+BCSTB={total:.4f} != BIOCHC={row['BIOCHC']:.4f}"
            )

        data_rows.append(row)

    if not data_rows:
        warnings.append("No data rows found — was biochar applied?")

    # Check: after first biochar application (NAPBIO>0), BIOCHC > 0
    applied_rows = [r for r in data_rows if r["NAPBIO"] > 0]
    if applied_rows:
        biochc_zero = [r for r in applied_rows if r["BIOCHC"] <= 0.0]
        if biochc_zero:
            errors.append(
                f"{len(biochc_zero)} rows with NAPBIO>0 but BIOCHC=0"
            )
    else:
        warnings.append("NAPBIO is always 0 — no biochar was applied this run")

    return errors, warnings


# ---------------------------------------------------------------------------
# Built-in self-test: synthetic BIOCHAR.OUT content
# ---------------------------------------------------------------------------
GOOD_CONTENT = """\
*BIOCHAR CARBON AND NITROGEN DYNAMICS
! Simulation run:    1
@YEAR DOY   DAS   BIOCHC    BCLBL    BCSTB   BIOCHN   BCMINC  DPH_L1  SORBNH4 NAPBIO    CUMBCC    CUMBCN
 2024  92     1     0.00     0.00     0.00     0.00     0.00     0.00     0.00      0      0.00      0.00
 2024  93     2  3000.00  150.00  2850.00   15.00     0.30     0.05     0.10      1   3000.00     15.00
 2024  94     3  2999.70  149.85  2849.85   14.99     0.30     0.05     0.10      1   3000.00     15.00
"""

BAD_CONTENT_MISMATCH = """\
*BIOCHAR CARBON AND NITROGEN DYNAMICS
@YEAR DOY   DAS   BIOCHC    BCLBL    BCSTB   BIOCHN   BCMINC  DPH_L1  SORBNH4 NAPBIO    CUMBCC    CUMBCN
 2024  92     1   100.00    10.00    80.00     0.50     0.01     0.00     0.00      1    100.00      0.50
"""


def self_test():
    ok = True

    def run(name, content, expect_errors, expect_warns):
        nonlocal ok
        fh = io.StringIO(content)
        lines = fh.readlines()
        tmp = "/tmp/_biochar_test_tmp.out"
        with open(tmp, "w") as f:
            f.writelines(lines)
        errs, warns = check_file(tmp)
        os.unlink(tmp)
        has_err = bool(errs)
        has_warn = bool(warns)
        passed = (has_err == expect_errors) and (has_warn == expect_warns)
        status = "PASS" if passed else "FAIL"
        if not passed:
            ok = False
            print(f"  {status} [{name}]  errors={errs}  warnings={warns}")
        else:
            print(f"  {status} [{name}]")

    print("Running BIOCHAR.OUT checker self-tests ...")
    run("good content",    GOOD_CONTENT,          expect_errors=False, expect_warns=False)
    run("mismatch BCLBL",  BAD_CONTENT_MISMATCH,  expect_errors=True,  expect_warns=False)
    print("All self-tests passed." if ok else "SELF-TEST FAILURES — see above.")
    return ok


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--self-test":
        success = self_test()
        sys.exit(0 if success else 1)

    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    path = sys.argv[1]
    errors, warnings = check_file(path)

    for w in warnings:
        print(f"WARNING: {w}")
    for e in errors:
        print(f"ERROR:   {e}")

    if errors:
        print(f"\nFAIL — {len(errors)} error(s) in {path}")
        sys.exit(1)
    else:
        print(f"PASS — {path} looks valid ({len(warnings)} warning(s))")
        sys.exit(0)
