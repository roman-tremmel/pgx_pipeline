// PyPGx (MIT). Star-allele + structural-variant caller.
//   NGS pipeline : WGS/WES BAM/CRAM  -> run-ngs-pipeline per gene
//   CHIP pipeline: array/imputed VCF -> run-chip-pipeline per gene
// Genes are iterated inside the container; results are merged to one TSV.

// Genes PyPGx has genotype/phenotype tables for and that we want in consensus.
def PYPGX_GENES = 'CYP2D6 CYP2C19 CYP2C9 CYP2B6 CYP3A5 DPYD TPMT NUDT15 UGT1A1 SLCO1B1 GSTM1 GSTT1 CYP2A6 NAT2'

process PYPGX_NGS {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(bai)
    val   assembly
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.pypgx.tsv"), emit: calls

    script:
    def cram = fasta ? "--reference ${fasta}" : ""
    """
    set -euo pipefail
    echo -e "gene\\tgenotype\\tphenotype\\ttool_detail" > ${meta.id}.pypgx.tsv
    for GENE in ${PYPGX_GENES}; do
        # per-gene NGS pipeline; tolerate genes with no model on this build
        if pypgx run-ngs-pipeline \$GENE ${meta.id}_\$GENE \\
                --variants <(echo) --assembly ${assembly} \\
                --bam ${bam} ${cram} 2> \$GENE.log; then
            if [ -f ${meta.id}_\$GENE/results.zip ]; then
                unzip -o -q ${meta.id}_\$GENE/results.zip -d ${meta.id}_\$GENE
                tail -n +2 ${meta.id}_\$GENE/data.tsv 2>/dev/null | \\
                    awk -v g=\$GENE 'BEGIN{OFS="\\t"}{print g,\$2,\$3,\$0}' \\
                    >> ${meta.id}.pypgx.tsv || true
            fi
        else
            echo "pypgx: skipped \$GENE (no model for ${assembly} or low coverage)" >&2
        fi
    done
    """

    stub:
    "echo -e 'gene\\tgenotype\\tphenotype\\ttool_detail\\nCYP2D6\\t*1/*4\\tIM\\tstub' > ${meta.id}.pypgx.tsv"
}

process PYPGX_CHIP {
    tag "$meta.id"

    input:
    tuple val(meta), path(vcf)
    val   assembly

    output:
    tuple val(meta), path("${meta.id}.pypgx.tsv"), emit: calls

    script:
    def genes_chip = 'CYP2C19 CYP2C9 CYP2B6 CYP3A5 DPYD TPMT NUDT15 UGT1A1 SLCO1B1 NAT2'
    """
    set -euo pipefail
    echo -e "gene\\tgenotype\\tphenotype\\ttool_detail" > ${meta.id}.pypgx.tsv
    for GENE in ${genes_chip}; do
        if pypgx run-chip-pipeline \$GENE ${meta.id}_\$GENE ${vcf} \\
                --assembly ${assembly} 2> \$GENE.log; then
            if [ -f ${meta.id}_\$GENE/results.zip ]; then
                unzip -o -q ${meta.id}_\$GENE/results.zip -d ${meta.id}_\$GENE
                tail -n +2 ${meta.id}_\$GENE/data.tsv 2>/dev/null | \\
                    awk -v g=\$GENE 'BEGIN{OFS="\\t"}{print g,\$2,\$3,\$0}' \\
                    >> ${meta.id}.pypgx.tsv || true
            fi
        fi
    done
    """

    stub:
    "echo -e 'gene\\tgenotype\\tphenotype\\ttool_detail\\nCYP2C19\\t*1/*2\\tIM\\tstub' > ${meta.id}.pypgx.tsv"
}
