// PharmCAT (MPL-2.0). Runs the VCF preprocessor then the full pipeline.
// Optional HLA (and other) outside calls are merged in. Emits the matcher TSV
// which the parser reads for per-gene star-allele calls.

process PHARMCAT {
    tag "$meta.id"

    input:
    tuple val(meta), path(vcf), path(outside)   // outside = NO_FILE sentinel when absent

    output:
    tuple val(meta), path("${meta.id}*.report.tsv"),  emit: calls
    tuple val(meta), path("${meta.id}*.report.html"), emit: report, optional: true
    tuple val(meta), path("${meta.id}*.match.json"),  emit: match,  optional: true

    script:
    def outside_arg = (outside.name != 'NO_FILE') ? "-po ${outside}" : ""
    """
    set -euo pipefail
    # pharmcat_pipeline wraps preprocess + named-allele matcher + phenotyper + reporter
    pharmcat_pipeline ${vcf} \\
        ${outside_arg} \\
        -reporterJson -matcherHtml \\
        -o . -bf ${meta.id}
    # produce a flat TSV of gene -> diplotype for the parser
    python3 - ${meta.id} > ${meta.id}.report.tsv <<'PY'
import sys, json, glob
sample = sys.argv[1]
print("gene\\tdiplotype\\tphenotype")
for f in glob.glob(f"{sample}*.report.json") + glob.glob(f"{sample}*.phenotype.json"):
    try:
        data = json.load(open(f))
    except Exception:
        continue
    genes = data.get("genes", {})
    # PharmCAT nests by data source (CPIC/DPWG); walk defensively
    def walk(obj):
        if isinstance(obj, dict):
            if "gene" in obj and ("diplotypes" in obj or "sourceDiplotypes" in obj):
                dips = obj.get("sourceDiplotypes") or obj.get("diplotypes") or []
                for d in dips:
                    name = d.get("label") or d.get("name") or "/".join(d.get("alleles",[]))
                    pheno = (d.get("phenotypes") or [""])[0]
                    print(f'{obj["gene"]}\\t{name}\\t{pheno}')
            for v in obj.values(): walk(v)
        elif isinstance(obj, list):
            for v in obj: walk(v)
    walk(data)
    break
PY
    """

    stub:
    "echo -e 'gene\\tdiplotype\\tphenotype\\nCYP2C19\\t*1/*2\\tIM' > ${meta.id}.report.tsv"
}
