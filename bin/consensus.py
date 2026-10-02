#!/usr/bin/env python3
"""Collapse harmonized per-tool calls into one consensus diplotype per gene.

Consensus rule (per gene, per sample):
  - group each tool's canonical_diplotype
  - the consensus is the diplotype supported by the most tools (plurality)
  - agreement level:
        FULL      : all tools that called the gene agree
        MAJORITY  : >50% of calling tools agree
        TIE/NONE  : no plurality winner (reported, flagged for review)
        SINGLE    : only one tool called the gene

Anything that is not FULL is also written to the discordance report with each
tool's call, so differences are visible and can feed harmonization-table curation.
"""
import argparse
import csv
from collections import defaultdict, Counter


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--genes", required=True, help="comma-separated genes to attempt consensus on")
    ap.add_argument("--consensus", required=True)
    ap.add_argument("--discordance", required=True)
    args = ap.parse_args()

    target_genes = {g.strip().upper() for g in args.genes.split(",") if g.strip()}

    # (sample, gene) -> {tool: canonical_diplotype}
    calls = defaultdict(dict)
    sample = ""
    rows = []
    with open(args.input) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            rows.append(r)
            sample = r.get("sample", sample)
            gene = r.get("gene", "").upper()
            tool = r.get("tool", "")
            dip = r.get("canonical_diplotype", "") or r.get("diplotype", "")
            if gene and tool and dip:
                # keep the first non-empty call per tool (tools may emit multiple solutions)
                calls[(sample, gene)].setdefault(tool, dip)

    cons_rows = []
    disc_rows = []
    for (smp, gene), tool_calls in sorted(calls.items()):
        # only form a consensus for requested genes; others are reported per-tool only
        attempt = gene in target_genes
        counts = Counter(tool_calls.values())
        n_tools = len(tool_calls)
        top_dip, top_n = counts.most_common(1)[0]
        tied = [d for d, c in counts.items() if c == top_n]

        if n_tools == 1:
            agreement = "SINGLE"
            winner = top_dip
        elif len(counts) == 1:
            agreement = "FULL"
            winner = top_dip
        elif len(tied) > 1:
            agreement = "TIE"
            winner = ""      # no winner; needs review
        elif top_n * 2 > n_tools:
            agreement = "MAJORITY"
            winner = top_dip
        else:
            agreement = "NONE"
            winner = ""

        supporting = ",".join(sorted(t for t, d in tool_calls.items() if d == winner)) if winner else ""

        if attempt:
            cons_rows.append({
                "sample": smp, "gene": gene, "consensus_diplotype": winner,
                "agreement": agreement, "n_tools": n_tools,
                "supporting_tools": supporting,
                "all_calls": "; ".join(f"{t}={d}" for t, d in sorted(tool_calls.items())),
            })

        if agreement != "FULL":
            for t, d in sorted(tool_calls.items()):
                disc_rows.append({"sample": smp, "gene": gene, "tool": t,
                                  "canonical_diplotype": d, "agreement": agreement})

    with open(args.consensus, "w", newline="") as fh:
        cols = ["sample", "gene", "consensus_diplotype", "agreement",
                "n_tools", "supporting_tools", "all_calls"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(cons_rows)

    with open(args.discordance, "w", newline="") as fh:
        cols = ["sample", "gene", "tool", "canonical_diplotype", "agreement"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(disc_rows)


if __name__ == "__main__":
    main()
