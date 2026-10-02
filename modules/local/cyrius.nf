// Cyrius (Illumina, Apache-2.0). CYP2D6 only, short-read WGS only.
// Consumes a manifest listing the BAM/CRAM; emits a per-sample TSV+JSON.

process CYRIUS {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(bai)
    val   genome                         // 19 | 37 | 38
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.cyrius.tsv"), emit: calls
    tuple val(meta), path("*.json"),                emit: json, optional: true

    script:
    def ref = fasta ? "--reference ${fasta}" : ""
    """
    set -euo pipefail
    echo "${bam}" > manifest.txt
    star_caller.py --manifest manifest.txt \\
        --genome ${genome} ${ref} \\
        --prefix ${meta.id} --outDir . --threads ${task.cpus}
    # normalize output name
    cp ${meta.id}.tsv ${meta.id}.cyrius.tsv
    """

    stub:
    "echo -e 'Sample\\tCYP2D6\\tFilter\\n${meta.id}\\t*1/*4\\tPASS' > ${meta.id}.cyrius.tsv"
}
