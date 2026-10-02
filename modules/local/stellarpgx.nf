// StellarPGx (MIT). Short-read high-coverage WGS. Normally its own Nextflow
// workflow; here we invoke the per-sample caller on a single BAM so it slots
// into this pipeline's channels. Emits one result file per sample.

def STELLAR_GENES = 'cyp2d6 cyp2c19 cyp2c9 cyp2b6 cyp3a5 cyp2a6 cyp3a4 gstm1 gstt1 nat2 slco1b1 nudt15 tpmt ugt1a1'

process STELLARPGX {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(bai)
    val   build                          // hg38 | b37 | hg19
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.stellarpgx.txt"), emit: calls

    script:
    """
    set -euo pipefail
    : > ${meta.id}.stellarpgx.txt
    for GENE in ${STELLAR_GENES}; do
        # StellarPGx exposes a core caller per gene; wrap and tag the output
        stellarpgx_call --gene \$GENE --build ${build} \\
            --bam ${bam} --ref ${fasta} --out ${meta.id}.\$GENE.out 2> \$GENE.log || \\
            echo "stellarpgx: \$GENE unavailable" >&2
        if [ -f ${meta.id}.\$GENE.out ]; then
            echo "## \$GENE" >> ${meta.id}.stellarpgx.txt
            cat ${meta.id}.\$GENE.out >> ${meta.id}.stellarpgx.txt
        fi
    done
    """

    stub:
    "printf '## CYP2D6\\nResult: *1/*4\\nActivity score: 1.0\\nMetaboliser status: IM\\n' > ${meta.id}.stellarpgx.txt"
}
