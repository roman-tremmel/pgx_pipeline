// PyPGx (MIT). Star-allele + structural-variant caller.
//   NGS pipeline : WGS/WES BAM/CRAM  -> run-ngs-pipeline per gene
//   CHIP pipeline: array/imputed VCF -> run-chip-pipeline per gene
// Every gene PyPGx implements is run; results are merged to one TSV.
// The gene list is enumerated from the installed PyPGx at runtime, so it tracks
// whatever the pinned container version supports.

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
    # all callable (target) genes implemented by this PyPGx build
    GENES=\$(python -c "import pypgx; print(' '.join(pypgx.list_genes(mode='target')))")
    echo -e "gene\\tgenotype\\tphenotype\\ttool_detail" > ${meta.id}.pypgx.tsv
    for GENE in \$GENES; do
        if pypgx run-ngs-pipeline \$GENE ${meta.id}_\$GENE \\
                --assembly ${assembly} \\
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
    """
    set -euo pipefail
    GENES=\$(python -c "import pypgx; print(' '.join(pypgx.list_genes(mode='target')))")
    echo -e "gene\\tgenotype\\tphenotype\\ttool_detail" > ${meta.id}.pypgx.tsv
    for GENE in \$GENES; do
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
