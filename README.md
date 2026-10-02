# pgx_pipeline

A Nextflow pipeline that runs several pharmacogenomics (PGx) star-allele callers
over the same samples, parses each tool's output into a common schema, harmonizes
star-allele nomenclature, and emits a **consensus call set** plus a **discordance
report**.

Different tools implement different genes and variants and spell the same allele
differently (e.g. a GSTM1 deletion as `del`, `*0`, or `-`). This pipeline makes
those calls comparable and records where the tools disagree.

> **Usage context:** this deployment is run for **academic, non-profit research**.
> That matters for Aldy's license (see [below](#aldy-licensing)).

---

## Contents

- [Tools](#tools)
- [Implementation status](#implementation-status)
- [Requirements](#requirements)
- [Quick start (tutorial)](#quick-start-tutorial)
- [The samplesheet](#the-samplesheet)
- [Command reference](#command-reference)
- [How harmonization works](#how-harmonization-works)
- [Consensus](#consensus)
- [Outputs](#outputs)

---

## Tools

| Tool | Scope | Data | Build | License | Default |
|------|-------|------|-------|---------|---------|
| [PharmCAT](https://pharmcat.org) 3.4.0 | CPIC genotype→phenotype→drug, star alleles from VCF | VCF (called or imputed) | 37/38 | MPL-2.0 | on |
| [PyPGx](https://pypgx.readthedocs.io) 0.26.0 | 87 genes, SV detection (del/dup/hybrid) | WGS/WES BAM, array VCF (chip pipeline) | 37/38 | MIT | on |
| [Cyrius](https://github.com/Illumina/Cyrius) | CYP2D6 only, highest WGS accuracy | WGS BAM | 19/37/38 | Apache-2.0 | on |
| [StellarPGx](https://github.com/SBIMB/StellarPGx) | ~20 genes incl. GSTM1/GSTT1 | WGS BAM | hg19/b37/hg38 | MIT | on |
| [OptiType](https://github.com/FRED-2/OptiType) | HLA-A/B (→ PharmCAT outside calls) | WGS/WES BAM | 37/38 | BSD-3 | on |
| [Aldy](https://github.com/0xTCG/aldy) 4.8.3 | 40+ genes, CN/phasing | WGS/WES BAM | hg19/hg38 | **Non-commercial (IU)** | on* |

### Aldy licensing

[Aldy](https://github.com/0xTCG/aldy) is licensed by Indiana University for
**non-commercial research only**. It **may not be used for commercial or
clinical/diagnostic purposes** without a separate license from IURTC.

Because this pipeline is used for **academic non-profit research**, Aldy is
enabled by default here (`run_aldy = true`, `aldy_license_ack = true` in
`nextflow.config`). The code keeps a safety gate: if `--run_aldy` is set without
`--aldy_license_ack`, the run aborts with the license notice.

**If you reuse this pipeline commercially or clinically**, set `--run_aldy false`
(or obtain a commercial license from Indiana University). All other tools
(PharmCAT, PyPGx, Cyrius, StellarPGx, OptiType) are usable commercially and
clinically under their respective OSI licenses.

### Data types

Short-read **WGS**, **WES**, and **genotyping arrays with imputation** are
supported; each caller runs only where it is valid:

- **WGS** → all callers (Cyrius and StellarPGx are WGS-only).
- **WES** → PyPGx (NGS), Aldy, OptiType, PharmCAT (called VCF).
- **array** → PyPGx (chip pipeline) and PharmCAT, from the imputed VCF.

Long-read-only callers are intentionally excluded.

**Every run calls every gene each tool implements** — PyPGx and Aldy enumerate
their gene lists from the installed container at runtime, StellarPGx runs its full
gene set, and Cyrius/OptiType cover their fixed targets (CYP2D6, HLA-A/B). The
`--consensus_genes` option only narrows the *consensus report*, not what the tools
call.

---

## Implementation status

All six tools are wired into the workflow: each has a Nextflow process, a parser
into the common schema, a pinned container image, and resource labels, and each
is routed by data type and genome build.

| Tool | Process module | Parser | Routed by datatype | Container command validated |
|------|----------------|--------|--------------------|-----------------------------|
| PharmCAT | `modules/local/pharmcat.nf` | `parse_pharmcat` | ✅ | ⚠️ not yet |
| PyPGx | `modules/local/pypgx.nf` (NGS + chip) | `parse_pypgx` | ✅ | ⚠️ not yet |
| Cyrius | `modules/local/cyrius.nf` | `parse_cyrius` | ✅ WGS only | ⚠️ not yet |
| StellarPGx | `modules/local/stellarpgx.nf` | `parse_stellarpgx` | ✅ WGS only | ⚠️ not yet |
| OptiType | `modules/local/optitype.nf` | `parse_optitype` | ✅ | ⚠️ not yet |
| Aldy | `modules/local/aldy.nf` | `parse_aldy` | ✅ | ⚠️ not yet |

**Tested:** the post-processing chain `parse_calls.py → harmonize.py →
consensus.py` is tested on synthetic multi-tool data.

**Not yet validated:** the per-tool container command lines have not been run
against the real CLIs. Expect to adjust: PyPGx's per-gene `run-ngs-pipeline` loop
and its `results.zip`/`data.tsv` layout; StellarPGx (normally its own Nextflow
workflow — the per-gene wrapper here is a placeholder); PharmCAT's report-JSON →
TSV extraction; and DeepVariant for PharmCAT input. Validate with the stub profile
first (below), then on one real sample per data type.

---

## Requirements

- [Nextflow](https://www.nextflow.io) ≥ 23.10 (needs Java 17+)
- [Singularity / Apptainer](https://apptainer.org) (the supported container engine;
  a `docker` profile is also provided)
- A reference FASTA matching `--genome_build` (required for CRAM input and for
  variant calling)

No tool needs a separate install — every caller runs from its container image,
pulled automatically on first use.

---

## Quick start (tutorial)

### 1. Get the pipeline

```bash
git clone https://github.com/roman-tremmel/pgx_pipeline.git
cd pgx_pipeline
```

### 2. Smoke-test the wiring (no data, no tool containers)

This runs every process with a `-stub` (a tiny fake command) so you can confirm
the workflow graph, channels, and post-processing all connect before touching real
data or pulling images:

```bash
nextflow run main.nf -profile test,singularity -stub
```

Look in `results_test/consensus/` for the stubbed consensus and discordance TSVs.

### 3. Point a cache at your Singularity image store (recommended)

```bash
export NXF_SINGULARITY_CACHEDIR=/path/to/shared/singularity_cache
```

### 4. Write a samplesheet

Copy the example and edit the paths:

```bash
cp assets/samplesheet.example.csv my_samples.csv
```

```csv
sample,datatype,bam,bai,vcf
HG001_wgs,wgs,/data/HG001.bam,/data/HG001.bam.bai,
NA12878_wes,wes,/data/NA12878.wes.bam,/data/NA12878.wes.bam.bai,
PATIENT_array,array,,,/data/PATIENT.imputed.vcf.gz
```

### 5. Run on real data

```bash
nextflow run main.nf \
    -profile singularity \
    --input my_samples.csv \
    --outdir results \
    --genome_build GRCh38 \
    --fasta /ref/GRCh38.fa
```

On an HPC scheduler add the `slurm` profile:

```bash
nextflow run main.nf -profile singularity,slurm \
    --input my_samples.csv --outdir results \
    --genome_build GRCh38 --fasta /ref/GRCh38.fa
```

### 6. Read the results

```bash
cat results/consensus/HG001_wgs.consensus.tsv        # one consensus diplotype per gene
cat results/consensus/HG001_wgs.discordance.tsv      # where tools disagreed
cat results/harmonize/HG001_wgs.unmapped.tsv         # allele names needing curation
```

### 7. Curate and re-run

Add a row to `assets/harmonization.tsv` for any allele in `*.unmapped.tsv`, then
re-run — Nextflow resumes with `-resume` and only the post-processing re-executes:

```bash
nextflow run main.nf -profile singularity -resume \
    --input my_samples.csv --outdir results \
    --genome_build GRCh38 --fasta /ref/GRCh38.fa
```

---

## The samplesheet

CSV with a header: `sample,datatype,bam,bai,vcf`

| column | required | notes |
|--------|----------|-------|
| `sample` | yes | unique sample id |
| `datatype` | yes | `wgs` \| `wes` \| `array` |
| `bam` | wgs/wes | path to BAM/CRAM |
| `bai` | optional | index; defaults to `<bam>.bai` if omitted |
| `vcf` | array (required); wgs/wes (optional) | for `array`, the imputed VCF. For wgs/wes, if given, PharmCAT uses it and variant calling is skipped |

---

## Command reference

Run the pipeline with `nextflow run main.nf [nextflow-options] [--pipeline-params]`.

### Built-in help

```bash
nextflow run main.nf --help
```

### Pipeline parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | — | **required.** Samplesheet CSV |
| `--outdir` | `results` | Output directory |
| `--genome_build` | `GRCh38` | `GRCh38` or `GRCh37` |
| `--fasta` | — | Reference FASTA (required for CRAM / variant calling) |
| `--fasta_fai` | `<fasta>.fai` | FASTA index |
| `--run_pypgx` | `true` | Toggle PyPGx |
| `--run_cyrius` | `true` | Toggle Cyrius (CYP2D6, WGS only) |
| `--run_stellarpgx` | `true` | Toggle StellarPGx (WGS only) |
| `--run_pharmcat` | `true` | Toggle PharmCAT |
| `--run_hla` | `true` | Toggle OptiType HLA typing |
| `--run_aldy` | `true` | Toggle Aldy (see licensing) |
| `--aldy_license_ack` | `true` | Must be true to run Aldy |
| `--consensus_genes` | `all` | Genes a consensus is computed for (`all` = every gene any tool reported). Each caller runs its full implemented gene set regardless |
| `--harmonization_table` | `assets/harmonization.tsv` | Curated lookup table |
| `--max_cpus` / `--max_memory` / `--max_time` | 16 / 128.GB / 48.h | Resource caps |

### Useful Nextflow options

| Option | Description |
|--------|-------------|
| `-profile singularity` | Use Singularity containers (add `,slurm` for HPC, `,test` for the test preset) |
| `-stub` | Run fake commands only — validate wiring without data or images |
| `-resume` | Reuse cached results; re-run only what changed |
| `-work-dir <dir>` | Scratch directory for intermediate work |
| `-with-report` / `-with-trace` | Extra run reports (also written to `pipeline_info/` by default) |

### Common recipes

```bash
# WGS cohort, GRCh38, all callers
nextflow run main.nf -profile singularity \
  --input cohort.csv --genome_build GRCh38 --fasta /ref/GRCh38.fa

# GRCh37 data
nextflow run main.nf -profile singularity \
  --input cohort.csv --genome_build GRCh37 --fasta /ref/hs37d5.fa

# Disable Aldy (e.g. for commercial/clinical reuse)
nextflow run main.nf -profile singularity \
  --input cohort.csv --fasta /ref/GRCh38.fa --run_aldy false

# Only PharmCAT + PyPGx
nextflow run main.nf -profile singularity \
  --input cohort.csv --fasta /ref/GRCh38.fa \
  --run_cyrius false --run_stellarpgx false --run_hla false --run_aldy false

# Restrict the CONSENSUS report to a gene subset (tools still call all genes)
nextflow run main.nf -profile singularity \
  --input cohort.csv --fasta /ref/GRCh38.fa \
  --consensus_genes "CYP2D6,CYP2C19,DPYD,TPMT"
```

---

## How harmonization works

Two layers (full detail in [docs/harmonization.md](docs/harmonization.md)):

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
`<sample>.unmapped.tsv` so a curator knows exactly what to add.

---

## Consensus

Per gene, per sample, the plurality canonical diplotype wins. Agreement is
labelled `FULL` / `MAJORITY` / `TIE` / `NONE` / `SINGLE`. Anything not `FULL`
is written to `<sample>.discordance.tsv`.

---

## Outputs

```
results/
├── pypgx/ aldy/ cyrius/ stellarpgx/ hla/ pharmcat/ variants/   # raw per-tool
├── parse_calls/                 # common-schema long tables
├── harmonize/  <sample>.harmonized.tsv  <sample>.unmapped.tsv
├── consensus/  <sample>.consensus.tsv   <sample>.discordance.tsv
└── pipeline_info/               # Nextflow timeline/report/trace/dag
```
