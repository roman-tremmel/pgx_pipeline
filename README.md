# pgx_pipeline

A Nextflow pipeline that runs several pharmacogenomics (PGx) star-allele callers
over the same samples, parses each tool's output into a common schema, harmonizes
star-allele nomenclature, and emits a **consensus call set** plus a **discordance
report**.

Different tools implement different genes and variants and spell the same allele
differently (e.g. a GSTM1 deletion as `del`, `*0`, or `-`). This pipeline makes
those calls comparable and records where the tools disagree.

## Tools

| Tool | Scope | Data | Build | License | Default |
|------|-------|------|-------|---------|---------|
| [PharmCAT](https://pharmcat.org) 3.4.0 | CPIC genotype→phenotype→drug, star alleles from VCF | VCF (called or imputed) | 37/38 | MPL-2.0 | on |
| [PyPGx](https://pypgx.readthedocs.io) 0.26.0 | 87 genes, SV detection (del/dup/hybrid) | WGS/WES BAM, array VCF (chip pipeline) | 37/38 | MIT | on |
| [Cyrius](https://github.com/Illumina/Cyrius) | CYP2D6 only, highest WGS accuracy | WGS BAM | 19/37/38 | Apache-2.0 | on |
| [StellarPGx](https://github.com/SBIMB/StellarPGx) | ~20 genes incl. GSTM1/GSTT1 | WGS BAM | hg19/b37/hg38 | MIT | on |
| [OptiType](https://github.com/FRED-2/OptiType) | HLA-A/B (→ PharmCAT outside calls) | WGS/WES BAM | 37/38 | BSD-3 | on |
| [Aldy](https://github.com/0xTCG/aldy) 4.8.3 | 40+ genes, CN/phasing | WGS/WES BAM | hg19/hg38 | **Non-commercial (IU)** | **off** |

### ⚠️ Aldy licensing

Aldy is licensed by Indiana University for **non-commercial research only**. It
**may not be used for commercial or clinical/diagnostic purposes** without a
separate license from IURTC. It is disabled by default and refuses to run unless
you assert eligibility:

```
--run_aldy --aldy_license_ack
```

All other tools are usable in commercial and clinical settings under their
respective OSI licenses.

### Data types

Short-read **WGS**, **WES**, and **genotyping arrays with imputation** are
supported; each caller runs only where it is valid:

- **WGS** → all callers (Cyrius and StellarPGx are WGS-only).
- **WES** → PyPGx (NGS), Aldy, OptiType, PharmCAT (called VCF).
- **array** → PyPGx (chip pipeline) and PharmCAT, from the imputed VCF.

Long-read-only callers are intentionally excluded.

## Usage

```bash
nextflow run main.nf -profile singularity \
    --input samplesheet.csv \
    --outdir results \
    --genome_build GRCh38 \
    --fasta /ref/GRCh38.fa
```

Smoke test (no data/containers, uses process stubs):

```bash
nextflow run main.nf -profile test,singularity -stub
```

### Samplesheet (`assets/samplesheet.example.csv`)

```csv
sample,datatype,bam,bai,vcf
HG001_wgs,wgs,/data/HG001.bam,/data/HG001.bam.bai,
NA12878_wes,wes,/data/NA12878.wes.bam,/data/NA12878.wes.bam.bai,
SAMPLE_array,array,,,/data/SAMPLE.imputed.vcf.gz
```

`datatype` is `wgs | wes | array`. For `wgs`/`wes` give a BAM/CRAM (`vcf`
optional — if present, PharmCAT uses it and variant calling is skipped). For
`array` give the imputed `vcf`.

### Key parameters

| Param | Default | Notes |
|-------|---------|-------|
| `--genome_build` | `GRCh38` | or `GRCh37` |
| `--fasta` | — | required for CRAM and for variant calling |
| `--run_pypgx/cyrius/stellarpgx/pharmcat/hla` | `true` | toggle callers |
| `--run_aldy` + `--aldy_license_ack` | `false` | see licensing above |
| `--consensus_genes` | CYP2D6,CYP2C19,… | genes a consensus is computed for |
| `--harmonization_table` | `assets/harmonization.tsv` | the curated lookup table |

## How harmonization works

Two layers (see [docs/harmonization.md](docs/harmonization.md)):

1. **Syntax rules** (regex, no curation): deletion spellings (`del`, `-`, `*0`,
   `0`, `null`) → a gene-specific canonical token (`*5` for CYP2D6, `*0` for
   GSTM1/GSTT1); duplication spellings (`*1x2`, `*1×2`) → `*1xN`.
2. **Curated lookup table** `assets/harmonization.tsv` — the file you grow over
   time. Columns:

   ```
   gene  tool  tool_allele  canonical_allele  relation  note  curator  date
   ```

   `tool` may be `*` (any tool). `relation` ∈ `synonym | subsumed_by | deletion |
   duplication | unsupported`. `subsumed_by` is for the case where one tool lacks
   an allele and calls a parent (e.g. `*1` when it means `*38`).

Alleles matching neither layer are passed through **and** listed in
`*.unmapped.tsv` so a curator knows exactly what to add to the table.

## Consensus

Per gene, per sample, the plurality canonical diplotype wins. Agreement is
labelled `FULL` / `MAJORITY` / `TIE` / `NONE` / `SINGLE`. Anything not `FULL`
is written to `*.discordance.tsv`.

## Outputs

```
results/
├── pypgx/ aldy/ cyrius/ stellarpgx/ hla/ pharmcat/ variants/   # raw per-tool
├── parse_calls/                 # common-schema long tables
├── harmonize/  <sample>.harmonized.tsv  <sample>.unmapped.tsv
├── consensus/  <sample>.consensus.tsv   <sample>.discordance.tsv
└── pipeline_info/               # Nextflow timeline/report/trace/dag
```

## Status

Scaffold. Workflow wiring, the common schema, harmonization + consensus logic and
their unit path are in place and tested on synthetic data. The per-tool container
command lines (marked in `modules/local/*.nf`) still need validating against each
tool's real CLI and output layout on a test dataset, and the container image tags
in `conf/modules.config` should be pinned to what your registry serves.
