#!/usr/bin/env python3
"""Harmonize per-tool star-allele nomenclature into a canonical form.

Two layers, applied in order:

  1. SYNTAX RULES (regex)   - mechanical fixes that need no curation:
       * deletion spellings:  del, -, *0, *5(CYP2D6), "deletion", "0"  -> gene default
       * duplication spellings: *1x2, *1X2, *1×2  -> *1xN  (copy number normalized)
       * whitespace / unicode x
  2. CURATED LOOKUP TABLE   - assets/harmonization.tsv, columns:
       gene  tool  tool_allele  canonical_allele  relation  note  curator  date
     relation in {synonym, subsumed_by, deletion, duplication, unsupported}
     tool may be '*' to match any tool. This is the file that grows over time.

Alleles that match neither layer are passed through unchanged AND written to the
unmapped report so a curator can add a row to the table.
"""
import argparse
import csv
import re
import sys
from collections import OrderedDict

LONG_COLUMNS = ["sample", "gene", "tool", "diplotype", "allele1", "allele2",
                "copy_number", "phenotype", "quality", "raw"]

# Gene-specific canonical token for a whole-gene deletion.
DELETION_CANONICAL = {
    "CYP2D6": "*5",       # PharmVar: *5 is the gene deletion
    "GSTM1": "*0",        # no PharmVar nomenclature; *0 = null/deletion
    "GSTT1": "*0",
    "UGT2B17": "*0",
    "_default": "*0",
}

DEL_TOKENS = {"del", "deletion", "-", "0", "*0", "null", "none"}


def normalize_syntax(gene, allele):
    """Layer 1: mechanical normalization. Returns (canonical_or_same, changed?)."""
    if allele is None:
        return "", False
    a = allele.strip()
    if a == "":
        return "", False

    low = a.lower()

    # --- deletions ---
    if low in DEL_TOKENS:
        return DELETION_CANONICAL.get(gene, DELETION_CANONICAL["_default"]), True

    # --- duplications: unify the multiplication sign and trailing count ---
    # *1x2, *1X2, *1×2, *1 x2  -> *1xN  (specific count preserved in copy_number col)
    dup = re.match(r"^(\*?\w+)\s*[x×X]\s*(\d+)$", a)
    if dup:
        return f"{dup.group(1)}xN", True

    # strip a stray leading 'rs'? leave rsIDs alone. normalize unicode star.
    a = a.replace("∗", "*")  # asterisk operator -> plain star
    return a, (a != allele.strip())


def load_table(path):
    """Return dict keyed by (gene, tool, tool_allele) -> (canonical, relation)."""
    tbl = {}
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            if not r.get("gene") or not r.get("tool_allele"):
                continue
            key = (r["gene"].strip(), r["tool"].strip(), r["tool_allele"].strip())
            tbl[key] = (r["canonical_allele"].strip(), r.get("relation", "synonym").strip())
    return tbl


def lookup(tbl, gene, tool, allele):
    """Exact (gene,tool,allele), then wildcard tool '*'. Returns (canonical, relation) or None."""
    for t in (tool, "*"):
        hit = tbl.get((gene, t, allele))
        if hit:
            return hit
    return None


def canonicalize_allele(tbl, gene, tool, allele, unmapped):
    """Apply syntax layer then table layer. Record unmapped non-trivial alleles."""
    syn, _ = normalize_syntax(gene, allele)
    if syn == "":
        return ""
    hit = lookup(tbl, gene, tool, syn)
    if hit:
        return hit[0]
    # also try the original (pre-syntax) spelling against the table
    hit = lookup(tbl, gene, tool, allele.strip())
    if hit:
        return hit[0]
    # unmapped: pass through, but flag anything that isn't already a clean *N
    if not re.match(r"^\*\w+(xN)?$", syn):
        unmapped.add((gene, tool, allele.strip(), syn))
    return syn


def canonical_diplotype(a1, a2):
    """Order-independent canonical string for comparison across tools."""
    alleles = [x for x in (a1, a2) if x]
    alleles.sort()
    return "/".join(alleles)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--unmapped", required=True)
    args = ap.parse_args()

    tbl = load_table(args.table)
    unmapped = set()
    out_rows = []

    for path in args.inputs:
        with open(path) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            for r in rd:
                gene = r.get("gene", "").strip()
                tool = r.get("tool", "").strip()
                if not gene:
                    continue
                c1 = canonicalize_allele(tbl, gene, tool, r.get("allele1", ""), unmapped)
                c2 = canonicalize_allele(tbl, gene, tool, r.get("allele2", ""), unmapped)
                r["canonical_allele1"] = c1
                r["canonical_allele2"] = c2
                r["canonical_diplotype"] = canonical_diplotype(c1, c2)
                out_rows.append(r)

    out_cols = LONG_COLUMNS + ["canonical_allele1", "canonical_allele2", "canonical_diplotype"]
    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

    with open(args.unmapped, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene", "tool", "tool_allele", "after_syntax_rules"])
        for g, t, a, s in sorted(unmapped):
            w.writerow([g, t, a, s])

    if unmapped:
        sys.stderr.write(f"harmonize: {len(unmapped)} allele spelling(s) need curation "
                         f"(see {args.unmapped})\n")


if __name__ == "__main__":
    main()
