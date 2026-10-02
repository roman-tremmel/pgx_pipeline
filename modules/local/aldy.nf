// Aldy (NON-COMMERCIAL / NON-DIAGNOSTIC, Indiana University license).
// Gated in main.nf behind --run_aldy --aldy_license_ack.
// WGS/WES BAM/CRAM. Runs the multi-gene genotyper and emits one combined file.

def ALDY_GENES = 'CYP2D6,CYP2C19,CYP2C9,CYP2B6,CYP3A5,CYP2A6,DPYD,TPMT,NUDT15,UGT1A1,SLCO1B1,NAT2'

process ALDY {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(bai)
    val   genome                         // hg19 | hg38
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.aldy"), emit: calls

    script:
    def ref = fasta ? "--reference ${fasta}" : ""
    """
    set -euo pipefail
    : > ${meta.id}.aldy
    for GENE in \$(echo ${ALDY_GENES} | tr ',' ' '); do
        aldy genotype -p wgs -g ${genome} ${ref} \\
            --gene \$GENE -o ${meta.id}.\$GENE.aldy ${bam} 2> \$GENE.log || \\
            echo "aldy: \$GENE failed/absent" >&2
        if [ -f ${meta.id}.\$GENE.aldy ]; then
            cat ${meta.id}.\$GENE.aldy >> ${meta.id}.aldy
        fi
    done
    """

    stub:
    "echo '#Solution 1: *1/*4' > ${meta.id}.aldy"
}
