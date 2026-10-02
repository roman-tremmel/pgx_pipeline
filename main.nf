#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
================================================================================
    pgx_pipeline
    Run several PGx star-allele callers per sample, parse each into a common
    schema, then harmonize nomenclature and emit a consensus call set.
================================================================================
*/

// --- caller modules ---
include { PYPGX_NGS   } from './modules/local/pypgx'
include { PYPGX_CHIP  } from './modules/local/pypgx'
include { ALDY        } from './modules/local/aldy'
include { CYRIUS      } from './modules/local/cyrius'
include { STELLARPGX  } from './modules/local/stellarpgx'
include { OPTITYPE    } from './modules/local/optitype'
include { DEEPVARIANT } from './modules/local/variantcalling'
include { PHARMCAT    } from './modules/local/pharmcat'

// --- post-processing ---
include { PARSE_CALLS } from './modules/local/parse'
include { HARMONIZE   } from './modules/local/harmonize'
include { CONSENSUS   } from './modules/local/consensus'

// ---------------------------------------------------------------------------
def helpMessage() {
    log.info """
    pgx_pipeline v${workflow.manifest.version}

    Usage:
      nextflow run main.nf -profile singularity \\
        --input samplesheet.csv --outdir results \\
        --genome_build GRCh38 --fasta genome.fa

    Samplesheet columns: sample,datatype,bam,bai,vcf
      datatype = wgs | wes | array
      wgs/wes  -> provide bam(+bai); vcf optional (skips variant calling for PharmCAT)
      array    -> provide an (imputed) vcf; bam columns left empty

    Key options:
      --genome_build   GRCh38 | GRCh37
      --run_aldy       enable Aldy (requires --aldy_license_ack, non-commercial research only)
      --run_cyrius / --run_stellarpgx / --run_pypgx / --run_pharmcat / --run_hla
    """.stripIndent()
}

// ---------------------------------------------------------------------------
workflow {

    if (params.help) { helpMessage(); return }
    if (!params.input) { error "You must supply --input samplesheet.csv" }

    // --- Aldy license gate (research item #4) -----------------------------
    if (params.run_aldy && !params.aldy_license_ack) {
        error """
        Aldy is licensed by Indiana University for NON-COMMERCIAL RESEARCH ONLY and may
        not be used for commercial or clinical/diagnostic purposes without a separate
        license from IURTC. To enable it you must assert eligibility:
            --run_aldy --aldy_license_ack
        """.stripIndent()
    }

    def build = params.genome_build
    if (!(build in ['GRCh38','GRCh37'])) { error "--genome_build must be GRCh38 or GRCh37" }
    def gmap  = params.genomes[build]

    // --- parse samplesheet into a typed channel ---------------------------
    ch_input = Channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [ id: row.sample, datatype: (row.datatype ?: 'wgs').toLowerCase() ]
            if (!(meta.datatype in ['wgs','wes','array']))
                error "sample ${meta.id}: datatype must be wgs|wes|array (got '${meta.datatype}')"
            [ meta, row ]
        }

    // alignment-based samples (wgs/wes) carry a BAM/CRAM
    ch_aln = ch_input
        .filter { meta, row -> meta.datatype in ['wgs','wes'] }
        .map    { meta, row ->
            if (!row.bam) error "sample ${meta.id}: datatype ${meta.datatype} requires a 'bam' file"
            tuple(meta, file(row.bam), row.bai ? file(row.bai) : file("${row.bam}.bai"))
        }

    ch_wgs = ch_aln.filter { meta, bam, bai -> meta.datatype == 'wgs' }

    // array samples carry an imputed VCF only
    ch_array_vcf = ch_input
        .filter { meta, row -> meta.datatype == 'array' }
        .map    { meta, row ->
            if (!row.vcf) error "sample ${meta.id}: datatype array requires a 'vcf' file"
            tuple(meta, file(row.vcf))
        }

    ch_fasta = params.fasta ? tuple(file(params.fasta),
                                    params.fasta_fai ? file(params.fasta_fai) : file("${params.fasta}.fai"))
                            : tuple([],[])

    // collects every parseable per-caller output: tuple(meta, tool, file)
    ch_raw = Channel.empty()

    // ---------------------- CALLERS ---------------------------------------
    if (params.run_pypgx) {
        // NGS pipeline for wgs/wes alignments; chip pipeline for array VCFs.
        // Each runs only for the samples in its (possibly empty) input channel.
        PYPGX_NGS  ( ch_aln,       gmap.pypgx_assembly, ch_fasta )
        PYPGX_CHIP ( ch_array_vcf, gmap.pypgx_assembly )
        ch_raw = ch_raw
            .mix( PYPGX_NGS.out.calls.map  { m,f -> tuple(m,'pypgx',f) } )
            .mix( PYPGX_CHIP.out.calls.map { m,f -> tuple(m,'pypgx',f) } )
    }

    if (params.run_aldy) {
        ALDY ( ch_aln, gmap.aldy_genome, ch_fasta )
        ch_raw = ch_raw.mix( ALDY.out.calls.map { m,f -> tuple(m,'aldy',f) } )
    }

    if (params.run_cyrius) {
        // CYP2D6 only, WGS only
        CYRIUS ( ch_wgs, gmap.cyrius_genome, ch_fasta )
        ch_raw = ch_raw.mix( CYRIUS.out.calls.map { m,f -> tuple(m,'cyrius',f) } )
    }

    if (params.run_stellarpgx) {
        STELLARPGX ( ch_wgs, gmap.stellarpgx_build, ch_fasta )
        ch_raw = ch_raw.mix( STELLARPGX.out.calls.map { m,f -> tuple(m,'stellarpgx',f) } )
    }

    if (params.run_hla) {
        OPTITYPE ( ch_aln )
        ch_raw = ch_raw.mix( OPTITYPE.out.calls.map { m,f -> tuple(m,'optitype',f) } )
    }

    // ---------------------- PHARMCAT --------------------------------------
    if (params.run_pharmcat) {
        def no_file = file("${projectDir}/assets/NO_FILE")

        // VCF supplied in the samplesheet (array samples, or wgs/wes rows with a vcf)
        ch_provided_vcf = ch_input
            .filter { meta, row -> row.vcf }
            .map    { meta, row -> tuple(meta, file(row.vcf)) }

        // alignment samples with NO vcf -> call variants. Each row is self-describing,
        // so this needs no join.
        ch_call = ch_input
            .filter { meta, row -> meta.datatype in ['wgs','wes'] && !row.vcf }
            .map    { meta, row ->
                tuple(meta, file(row.bam), row.bai ? file(row.bai) : file("${row.bam}.bai")) }

        DEEPVARIANT ( ch_call, ch_fasta )

        ch_pharmcat_vcf = DEEPVARIANT.out.vcf.mix( ch_provided_vcf )

        // attach optional HLA outside-call file by sample id; absent -> sentinel
        ch_outside = ( params.run_hla
            ? OPTITYPE.out.outside_call.map { meta, f -> tuple(meta.id, f) }
            : Channel.empty() )

        ch_pharmcat_in = ch_pharmcat_vcf
            .map  { meta, vcf -> tuple(meta.id, meta, vcf) }
            .join ( ch_outside, remainder: true )
            .map  { id, meta, vcf, outside -> tuple(meta, vcf, outside ?: no_file) }

        PHARMCAT ( ch_pharmcat_in )
        ch_raw = ch_raw.mix( PHARMCAT.out.calls.map { m,f -> tuple(m,'pharmcat',f) } )
    }

    // ---------------------- PARSE / HARMONIZE / CONSENSUS -----------------
    PARSE_CALLS ( ch_raw )

    // gather all per-sample long tables
    ch_long = PARSE_CALLS.out.tsv
        .map { meta, f -> tuple(meta.id, meta, f) }
        .groupTuple()
        .map { id, metas, files -> tuple(metas[0], files) }

    HARMONIZE ( ch_long, file(params.harmonization_table) )
    CONSENSUS ( HARMONIZE.out.harmonized, params.consensus_genes )
}

workflow.onComplete {
    log.info ( workflow.success
        ? "\npgx_pipeline complete. Results in ${params.outdir}\n"
        : "\npgx_pipeline finished with errors.\n" )
}
