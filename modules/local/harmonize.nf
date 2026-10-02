// Apply the harmonization layer (syntax rules + curated lookup table) to the
// combined per-sample long table. Adds canonical_allele1/2 and canonical_diplotype
// columns and flags unmapped alleles for curation.

process HARMONIZE {
    tag "$meta.id"

    input:
    tuple val(meta), path(long_tsvs)
    path  table

    output:
    tuple val(meta), path("${meta.id}.harmonized.tsv"),  emit: harmonized
    tuple val(meta), path("${meta.id}.unmapped.tsv"),    emit: unmapped

    script:
    """
    harmonize.py --table ${table} \\
        --inputs ${long_tsvs} \\
        --output ${meta.id}.harmonized.tsv \\
        --unmapped ${meta.id}.unmapped.tsv
    """

    stub:
    """
    echo -e 'sample\\tgene\\ttool\\tcanonical_diplotype' > ${meta.id}.harmonized.tsv
    echo -e 'gene\\ttool\\tallele' > ${meta.id}.unmapped.tsv
    """
}
