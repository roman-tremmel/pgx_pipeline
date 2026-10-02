#!/usr/bin/env python3
"""Parse a single caller's native output into the common long schema.

Common schema (TSV):
    sample  gene  tool  diplotype  allele1  allele2  copy_number  phenotype  quality  raw

Each tool has its own parser function. They are intentionally defensive: formats
drift between versions, so unknown lines are skipped rather than fatal, and the
original line is preserved in the `raw` column for debugging / curation.
"""
import argparse
import csv
import json
import os
import re
import sys

COLUMNS = ["sample", "gene", "tool", "diplotype", "allele1", "allele2",
           "copy_number", "phenotype", "quality", "raw"]


def split_diplotype(dip):
    """'*1/*4' -> ('*1', '*4'); tolerate missing/!=2 alleles."""
    if not dip:
        return "", ""
    parts = re.split(r"[/|]", dip.strip())
    if len(parts) == 2:
        return parts[0].strip(), parts[1].strip()
    if len(parts) == 1:
        return parts[0].strip(), ""
    return parts[0].strip(), parts[-1].strip()


def row(sample, tool, gene, dip, cn="", pheno="", qual="", raw=""):
    a1, a2 = split_diplotype(dip)
    return {
        "sample": sample, "gene": gene, "tool": tool, "diplotype": dip,
        "allele1": a1, "allele2": a2, "copy_number": cn,
        "phenotype": pheno, "quality": qual, "raw": raw.replace("\t", " ").strip(),
    }


# --------------------------------------------------------------------------
def parse_pypgx(path, sample, tool):
    out = []
    with open(path) as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd, None)
        for r in rd:
            if len(r) < 2:
                continue
            gene, dip = r[0], r[1]
            pheno = r[2] if len(r) > 2 else ""
            out.append(row(sample, tool, gene, dip, pheno=pheno, raw="\t".join(r)))
    return out


def parse_cyrius(path, sample, tool):
    out = []
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            dip = r.get("CYP2D6") or r.get("Genotype") or ""
            filt = r.get("Filter", "")
            # Cyrius may report e.g. "*1/*4" or copy-number style "*1x2/*4"
            out.append(row(sample, tool, "CYP2D6", dip, qual=filt, raw=json.dumps(r)))
    return out


def parse_aldy(path, sample, tool):
    """Aldy .aldy files start solutions with '#Solution ...: *a/*b' and carry a
    gene name per block. We scan for solution lines and the gene header."""
    out = []
    gene = ""
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            m_gene = re.match(r"#\s*Gene:\s*(\S+)", line)
            if m_gene:
                gene = m_gene.group(1)
                continue
            m = re.match(r"#\s*Solution\s*\d+:\s*(.+)", line)
            if m:
                dip = m.group(1).strip()
                out.append(row(sample, tool, gene or "UNKNOWN", dip, raw=line))
    return out


def parse_stellarpgx(path, sample, tool):
    """Blocks headed by '## GENE' with lines like 'Result: *1/*4' and
    'Activity score: ...' / 'Metaboliser status: ...'."""
    out = []
    gene, dip, pheno = "", "", ""
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("##"):
                if gene and dip:
                    out.append(row(sample, tool, gene.upper(), dip, pheno=pheno))
                gene, dip, pheno = line.lstrip("#").strip(), "", ""
            elif line.lower().startswith("result"):
                dip = line.split(":", 1)[1].strip()
            elif "metaboliser" in line.lower() or "metabolizer" in line.lower():
                pheno = line.split(":", 1)[1].strip() if ":" in line else line
        if gene and dip:
            out.append(row(sample, tool, gene.upper(), dip, pheno=pheno))
    return out


def parse_pharmcat(path, sample, tool):
    out = []
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            gene = r.get("gene", "")
            dip = r.get("diplotype") or r.get("Diplotype") or ""
            pheno = r.get("phenotype", "")
            if gene:
                out.append(row(sample, tool, gene, dip, pheno=pheno, raw="\t".join(r.values())))
    return out


def parse_optitype(path, sample, tool):
    """Emit HLA-A and HLA-B as two 'genes' in the long schema."""
    out = []
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            for gene, cols in [("HLA-A", ["A1", "A2"]), ("HLA-B", ["B1", "B2"])]:
                a = [r.get(c, "") for c in cols if r.get(c)]
                if len(a) == 2:
                    out.append(row(sample, tool, gene, "/".join(a), raw=json.dumps(r)))
            break  # OptiType reports a single best row
    return out


PARSERS = {
    "pypgx": parse_pypgx,
    "cyrius": parse_cyrius,
    "aldy": parse_aldy,
    "stellarpgx": parse_stellarpgx,
    "pharmcat": parse_pharmcat,
    "optitype": parse_optitype,
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--tool", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    parser = PARSERS.get(args.tool)
    if parser is None:
        sys.exit(f"parse_calls: unknown tool '{args.tool}'")

    rows = []
    if os.path.getsize(args.input) > 0:
        try:
            rows = parser(args.input, args.sample, args.tool)
        except Exception as e:  # never fail the whole pipeline on one parse
            sys.stderr.write(f"parse_calls: {args.tool} parse error: {e}\n")

    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


if __name__ == "__main__":
    main()
