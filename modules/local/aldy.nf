// Aldy (NON-COMMERCIAL / NON-DIAGNOSTIC, Indiana University license).
// Enabled for this academic non-profit deployment; gated in main.nf behind
// --run_aldy --aldy_license_ack.
// WGS/WES BAM/CRAM. Every gene Aldy ships a definition for is genotyped; the
// gene list is read from the installed Aldy's bundled resources at runtime.

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
    # enumerate every gene Aldy implements (one YAML definition per gene)
    GENES=\$(python - <<'PY'
import aldy, os, glob
d = os.path.join(os.path.dirname(aldy.__file__), "resources", "genes")
print(" ".join(sorted(os.path.splitext(os.path.basename(f))[0]
                      for f in glob.glob(os.path.join(d, "*.yml")))))
PY
)
    : > ${meta.id}.aldy
    for GENE in \$GENES; do
        if aldy genotype -p wgs -g ${genome} ${ref} \\
                --gene \$GENE -o ${meta.id}.\$GENE.aldy ${bam} 2> \$GENE.log; then
            # tag each block with its gene so the parser can attribute solutions
            echo "# Gene: \$GENE" >> ${meta.id}.aldy
            cat ${meta.id}.\$GENE.aldy >> ${meta.id}.aldy
        else
            echo "aldy: \$GENE failed/absent" >&2
        fi
    done
    """

    stub:
    "printf '# Gene: CYP2D6\\n#Solution 1: *1/*4\\n' > ${meta.id}.aldy"
}
