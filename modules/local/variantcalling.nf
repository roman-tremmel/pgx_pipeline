// Variant calling for PharmCAT input when no VCF was supplied.
// DeepVariant over the whole BAM (or restrict to a PGx BED for speed).
// WGS uses model_type WGS; WES uses WES.

process DEEPVARIANT {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.vcf.gz"), emit: vcf

    when:
    fasta   // DeepVariant needs a reference

    script:
    def model = meta.datatype == 'wes' ? 'WES' : 'WGS'
    """
    set -euo pipefail
    run_deepvariant \\
        --model_type=${model} \\
        --ref=${fasta} \\
        --reads=${bam} \\
        --output_vcf=${meta.id}.vcf.gz \\
        --num_shards=${task.cpus}
    """

    stub:
    "echo '' | gzip > ${meta.id}.vcf.gz"
}
