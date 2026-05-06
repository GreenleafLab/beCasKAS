# beCasKAS

Reproducible Nextflow pipeline implementing the beCasKAS computational analysis from
[Wang T. et al., bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.09.25.678665v1).
Detects non-coding off-target base editor unwinding and editing in primary cells.

## Quick start

### Requirements
- [Nextflow](https://www.nextflow.io/) >= 24.04.0
- Pipeline tools on the executor's `PATH` (see [Tool dependencies](#tool-dependencies) -- list pending)
- ~50 GB scratch space + reference data (~30 GB for hg38 + GATK resources)

### Run

```bash
nextflow run GreenleafLab/beCasKAS -r v0.2.0 \
    --samplesheet samples.csv \
    --reference_dir /path/to/hg38_refs \
    --outdir results/
```

The pipeline does not download references -- you stage them once and pass the directory.

For SLURM clusters, use `-profile slurm` and submit a thin driver job; see [submit_example.sbatch](submit_example.sbatch) for a template.

## Container support

The pipeline ships an Apptainer recipe under [nextflow/containers/becaskas.def](nextflow/containers/becaskas.def). Build the SIF locally from the repository root:

```bash
apptainer build --fakeroot becaskas.sif nextflow/containers/becaskas.def
```

`--fakeroot` enables unprivileged builds and works on most modern Linux distributions (requires user namespace support; available since kernel 3.8). HPC sites that disable user namespaces will require the admin to build the SIF and distribute it.

Run the pipeline with a `process.container` pointing at the local SIF:

```bash
nextflow run main.nf \
    -c <(echo "process.container = '$PWD/becaskas.sif'") \
    --samplesheet samples.csv --reference_dir /path/to/refs --outdir results/
```

Apptainer is the default container runtime (configured in [nextflow.config](nextflow.config)) and works for both direct `nextflow run` invocations and `-profile slurm` (HPC). Docker and Singularity profiles are also pre-wired but require a separate image; registry hosting (publishing a prebuilt image) is deferred to a future release.

## Required references

Stage the following layout under `--reference_dir`:

```
${reference_dir}/
  hg38.fasta                              # UCSC-style hg38 primary FASTA
  hg38.fasta.fai                          # samtools faidx
  hg38.dict                               # samtools dict (sibling of .fasta, GATK convention)
  hg38.genome                             # chromosome sizes (UCSC chrom.sizes format)
  hg38.1.bt2                              # bowtie2 forward index (4 files)
  hg38.2.bt2
  hg38.3.bt2
  hg38.4.bt2
  hg38.rev.1.bt2                          # bowtie2 reverse index (2 files)
  hg38.rev.2.bt2
  af-only-gnomad.hg38.vcf.gz              # GATK germline resource
  af-only-gnomad.hg38.vcf.gz.tbi
  1000g_pon.hg38.vcf.gz                   # GATK panel of normals
  1000g_pon.hg38.vcf.gz.tbi
```

All reference files live flat in `${reference_dir}/` (no subdirectories). The pipeline expects exactly these filenames. Suggested staging steps:
1. Download UCSC `hg38.fa.gz` from `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz`, gunzip, rename to `hg38.fasta`.
2. `samtools faidx hg38.fasta` and `samtools dict hg38.fasta -o hg38.dict`.
3. `bowtie2-build hg38.fasta hg38` (~30 to 60 minutes); produces the six `hg38.*.bt2` files in the current directory.
4. `samtools view -H hg38.fasta` to derive `hg38.genome` (chrom.sizes), or download from UCSC.
5. Download GATK somatic resources from `https://storage.googleapis.com/gatk-best-practices/somatic-hg38/`.

## Samplesheet

CSV with columns: `sample_id`, `fastq_1`, `fastq_2`, `condition`, `replicate`, `control_id`. See [samplesheet_example.csv](samplesheet_example.csv).

| Column | Required | Description |
|---|---|---|
| `sample_id` | yes | Unique sample identifier |
| `fastq_1` | yes | Path/URL to R1 FASTQ |
| `fastq_2` | no | Path/URL to R2 FASTQ (omit for single-end) |
| `condition` | yes | One of `input`, `control`, `treatment` |
| `replicate` | yes | Replicate number within condition |
| `control_id` | conditional | sample_id of MACS2 input control (required for control/treatment rows) |

## Tool dependencies

The canonical tool/version list lives in:
- [nextflow/containers/env.yaml](nextflow/containers/env.yaml) -- conda environment spec (bioinformatics tools: BWA, samtools, bcftools, bedtools, MACS2, trim_galore, FastQC, deepTools, picard, gencore, GATK4, UCSC bedGraphToBigWig).
- [nextflow/containers/install-r-packages.R](nextflow/containers/install-r-packages.R) -- R/Bioconductor packages (DiffBind, DESeq2, edgeR, GenomicRanges, rtracklayer + supporting CRAN packages).
- KAS-Analyzer 1.0 is cloned from `https://github.com/Ruitulyu/KAS-Analyzer.git` (see [becaskas.def](nextflow/containers/becaskas.def) `%post` section).

These two files are the source of truth: `becaskas.def` consumes them at build time, so any pin update applies to both the container build and any local conda environment built from the same `env.yaml`.

## Pipeline structure

| Stage | Process | Source file | Source script | Description |
|---|---|---|---|---|
| 1 | `TRIM_ALIGN` | `trim_align.nf` | 1_trim_align.sh / 2_align_helper.sh | trim_galore + KAS-Analyzer KAS-seq (Bowtie2, KAS-Analyzer's default aligner) |
| 2 | `PEAK_CALL` | `analyze.nf` | 3_analyze.sh | MACS2 sharp via KAS-Analyzer (treatment vs input) |
| 3 | `KAS_STATISTICS` | `analyze.nf` | 3_analyze.sh | Cohort RPM statistics -> ratios.txt |
| 4 | `NORMALIZE_BIGWIG` | `analyze.nf` | 3_analyze.sh | RPM-ratio normalized bedGraph + UCSC bigWig |
| 5 | `QC_METRICS` | `analyze.nf` | 3_analyze.sh | FRiP + multiBigwigSummary correlation |
| 6 | `DEDUP_CONSENSUS` | `picard_gencore_dedup.nf` | 4_picard-gencore.sh | picard MarkDuplicates + gencore consensus dedup |
| 7 | `SPLIT_BAM_PSEUDOREPS` | `split_bam_pseudoreps.nf` | (new) | (conditional) deterministic BAM partition for n=1 conditions |
| 8 | `DIFFERENTIAL_PEAKS` | `differential_peaks.nf` | 5E_differential_peaks.Rmd | DiffBind (DESeq2 backend) treatment vs control contrast |
| 9 | `CALL_VARIANTS` | `mutect2_call_variants.nf` | 6_gatk-mutect2.sh | Mutect2 tumor-only + germline + PoN, restricted to differential peaks |

## Replicate handling

| Samples per condition | Behavior |
|---|---|
| n >= 2 | Standard DiffBind path with biological replicates (publication-grade) |
| n = 1 | Auto-split BAM into 3 deterministic pseudo-replicates (`--n_pseudo_replicates`); fold-change estimates remain valid; p-values are optimistic -- **not for biological inference** |
| n = 0 (non-input) | Pipeline fails fast with a clear error |

## Parameters

Run `nextflow run <pipeline> --help` for the full list with descriptions.

Key parameters and their defaults (matching the original preprint code):

| Param | Default | Role |
|---|---|---|
| `--clip_r2` | 12 | Bases clipped from 5' of R2 (adaptase library prep) |
| `--macs2_mode` | sharp | MACS2 peak-calling mode |
| `--fdr` | 0.05 | DiffBind FDR cutoff |
| `--lfc` | 1 | DiffBind log2 fold-change cutoff |
| `--min_overlap` | 2 | DiffBind consensus minOverlap |
| `--n_pseudo_replicates` | 3 | Pseudo-rep count when n=1 |
| `--seed` | 42 | Global random seed |

## Output

```
results/
  01_trim_align/${sample_id}/             # sorted.bam, .bed, .bg, fastqc reports
  02_peaks/${sample_id}/                   # MACS2 peak BED + xls
  03_bigwigs/${sample_id}/                 # *.nor.bigWig, *.nor.bg
  04_qc/                                   # FRiP report, correlation PDF/tab
  05_dedup/${sample_id}/                   # *_picard-gencore-dedup.bam
    pseudoreps/                            # (only if n=1) pseudo-rep BAMs
  06_differential/                         # differential_peaks.bed, full TSV, volcano, dba_object.rds, method.txt
  07_variants/${sample_id}/                # *.vcf, *_no-alt-chr.vcf, *_CT.vcf, *_AG.vcf
  pipeline_info/                           # execution_report, timeline, trace, DAG SVG
```

## Reproducibility

- Random seeds set explicitly in R (`set.seed`).
- References are user-staged at `--reference_dir`; provenance of the staged files is the user's responsibility.
- Tool versions captured per-run at workflow start, written to `.nextflow.log` (look for the `Tool versions for this run:` block).
- Two runs with identical inputs and the same staged references produce bit-identical outputs (subject to bowtie2 multi-thread caveat -- for strict bit-identity, pin `task.cpus = 1` for the trim_align process).
- Container-based execution is supported via a locally-built Apptainer SIF (see [Container support](#container-support)). Byte-pinned reproducibility (registry-hosted image with sha256 digest pin) is deferred to a future release; v0.1.0 reproducibility is at the "same recipe, same versions" level rather than "same digest".

**Note:** outputs from this pipeline may differ trivially from the original preprint figures, because random seeds were not set in the original code. The original preprint figures are reproducible only via the legacy repository at commit `d811c3a` (https://github.com/<legacy-org>/beCasKAS).

## Citation

See [CITATION.cff](CITATION.cff). If you use this pipeline, please cite the preprint:

> Wang T., Jessa S., Marinov G.K., Klemm S., Kundaje A., Greenleaf W.J. (2025).
> Sensitive, direct detection of non-coding off-target base editor unwinding and editing in primary cells.
> bioRxiv. doi:10.1101/2025.09.25.678665

## License

MIT -- see [LICENSE](LICENSE).
