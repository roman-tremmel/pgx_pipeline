// OptiType (BSD-3). HLA-A/B/C class-I typing from WGS/WES.
// We extract reads over the MHC region, then type. Emits the raw result TSV
// and a PharmCAT 'outside call' TSV (HLA-A, HLA-B) for the reporter.

process OPTITYPE {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.hla.tsv"),          emit: calls
    tuple val(meta), path("${meta.id}.outside_calls.tsv"), emit: outside_call

    script:
    // MHC coordinates differ by build; default to GRCh38 chr6 region here.
    def region = params.genome_build == 'GRCh38' ? 'chr6:29000000-34000000' : '6:29000000-34000000'
    """
    set -euo pipefail
    samtools view -b ${bam} ${region} > mhc.bam
    samtools fastq -1 r1.fq -2 r2.fq -0 /dev/null -s /dev/null -n mhc.bam
    OptiTypePipeline.py -i r1.fq r2.fq --dna -v -o optitype_out
    RES=\$(find optitype_out -name '*_result.tsv' | head -n1)
    cp \$RES ${meta.id}.hla.tsv

    # build PharmCAT outside-call lines from the A/B columns
    python3 - "\$RES" ${meta.id} > ${meta.id}.outside_calls.tsv <<'PY'
import sys, csv
res, sample = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(res), delimiter='\t'))
if rows:
    r = rows[0]
    def fmt(x):   # "A*01:01" -> "*01:01"
        return '*' + x.split('*',1)[1] if x and '*' in x else None
    for gene, cols in [('HLA-A',['A1','A2']), ('HLA-B',['B1','B2'])]:
        alleles = [fmt(r.get(c,'')) for c in cols if fmt(r.get(c,''))]
        if len(alleles) == 2:
            print(f"{gene}\t{alleles[0]}/{alleles[1]}")
PY
    """

    stub:
    """
    echo -e 'A1\\tA2\\tB1\\tB2' > ${meta.id}.hla.tsv
    echo -e 'HLA-B\\t*57:01/*07:02' > ${meta.id}.outside_calls.tsv
    """
}
