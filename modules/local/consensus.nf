// Collapse harmonized per-tool calls into one consensus diplotype per gene,
// with an agreement level and a discordance report.

process CONSENSUS {
    tag "$meta.id"

    input:
    tuple val(meta), path(harmonized)
    val   consensus_genes

    output:
    tuple val(meta), path("${meta.id}.consensus.tsv"),   emit: consensus
    tuple val(meta), path("${meta.id}.discordance.tsv"), emit: discordance

    script:
    """
    consensus.py --input ${harmonized} \\
        --genes "${consensus_genes}" \\
        --consensus ${meta.id}.consensus.tsv \\
        --discordance ${meta.id}.discordance.tsv
    """

    stub:
    """
    echo -e 'sample\\tgene\\tconsensus_diplotype\\tagreement\\tn_tools\\tsupporting_tools' > ${meta.id}.consensus.tsv
    echo -e 'sample\\tgene\\ttool\\tcanonical_diplotype' > ${meta.id}.discordance.tsv
    """
}
