// Parse each caller's native output into the common long schema:
//   sample, gene, tool, diplotype, allele1, allele2, copy_number, phenotype, quality, raw
// One invocation per (sample, tool) raw file.

process PARSE_CALLS {
    tag "${meta.id}:${tool}"

    input:
    tuple val(meta), val(tool), path(raw)

    output:
    tuple val(meta), path("${meta.id}.${tool}.long.tsv"), emit: tsv

    script:
    """
    parse_calls.py --sample ${meta.id} --tool ${tool} \\
        --input ${raw} --output ${meta.id}.${tool}.long.tsv
    """

    stub:
    "echo -e 'sample\\tgene\\ttool\\tdiplotype\\tallele1\\tallele2\\tcopy_number\\tphenotype\\tquality\\traw' > ${meta.id}.${tool}.long.tsv"
}
