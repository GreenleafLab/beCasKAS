# beCasKAS

Reproducible Nextflow pipeline implementing the beCasKAS computational analysis from
[Wang T. et al., bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.09.25.678665v1).
Detects non-coding off-target base editor unwinding and editing in primary cells.

## Quick start

### Requirements
- [Nextflow](https://www.nextflow.io/) >= 24.04.0
- [Apptainer](https://apptainer.org/docs/admin/main/installation.html) >= 1.2.0
- ~50 GB scratch space + reference data (~30 GB for hg38 + GATK resources)

### Run

```bash
nextflow run GreenleafLab/beCasKAS \
    --samplesheet samples.csv \
    --refdir /path/to/refs \
    --outdir results/
```

The pipeline does not download references per run. They must be staged in a directory and passed to `--refdir`.

For SLURM clusters, use `-profile slurm` and submit a thin driver job; see [submit_example.sbatch](submit_example.sbatch) for a template.

## Container support

The pipeline ships an Apptainer recipe under [env/becaskas.def](env/becaskas.def). Build the SIF locally from the repository root:

```bash
apptainer build --fakeroot becaskas.sif env/becaskas.def
```

- `--fakeroot` enables unprivileged builds and works on most modern Linux distributions (requires user namespace support). HPCs that disable user namespaces will require the admin to build the SIF and distribute it.
- HPCs may have login/interactive download and memory limits. It is typically advised to run a `sbatch` SLURM job to build the container (See [container_build_example.sbatch](container_build_example.sbatch))

By default the pipeline expects `becaskas.sif` at the repository root. To use a SIF in a different location, override the default:

```bash
nextflow run main.nf --container /abs/path/to/becaskas.sif ...
```

Apptainer is the default container runtime (configured in [nextflow.config](nextflow.config)) and works for both direct `nextflow run` invocations and `-profile slurm` (HPC). A Docker profile is pre-configured, but requires separate installation.

## Required references

Stage the following reference layout under `--refdir`:

```
hg38.fasta                              # UCSC hg38 fasta
hg38.fasta.fai                          # samtools faidx
hg38.dict                               # samtools dict
hg38.genome                             # chromosome sizes (UCSC chrom.sizes)
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

All reference files live flat in `${refdir}/`. The pipeline expects exactly these filenames. Suggested staging steps:
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
| `control_id` | conditional | sample_id of MACS2 input control |

## Tool dependencies

Tool dependencies and pinned versions live in:
- [env/env.yaml](env/env.yaml) - conda environment spec (bioinformatics tools: BWA, samtools, bcftools, bedtools, MACS2, trim_galore, FastQC, deepTools, picard, gencore, GATK4, UCSC bedGraphToBigWig).
- KAS-Analyzer 2.0.0 is cloned from `https://github.com/Ruitulyu/KAS-Analyzer.git`

## Pipeline structure

| Stage | Process | Source file | Description |
|---|---|---|---|
| 1 | `TRIM_ALIGN` | `trim_align.nf` | trim_galore + KAS-Analyzer KAS-seq (Bowtie2, KAS-Analyzer's default aligner) |
| 2 | `PEAK_CALL` | `analyze.nf` | MACS2 sharp via KAS-Analyzer (treatment vs input) |
| 3 | `KAS_STATISTICS` | `analyze.nf` | Cohort RPM statistics -> ratios.txt |
| 4 | `NORMALIZE_BIGWIG` | `analyze.nf` | RPM-ratio normalized bedGraph + UCSC bigWig |
| 5 | `QC_METRICS` | `analyze.nf` | FRiP + multiBigwigSummary correlation |
| 6 | `DEDUP_CONSENSUS` | `picard_gencore_dedup.nf` | picard MarkDuplicates + gencore consensus dedup |
| 7 | `SPLIT_BAM_PSEUDOREPS` | `split_bam_pseudoreps.nf` | (conditional) deterministic BAM partition for n=1 conditions |
| 8 | `DIFFERENTIAL_PEAKS` | `differential_peaks.nf` | DiffBind (DESeq2 backend) treatment vs control contrast |
| 9 | `CALL_VARIANTS` | `mutect2_call_variants.nf` | Mutect2 tumor-only + germline + PoN, restricted to differential peaks |

## Replicate handling

| Samples per condition | Behavior |
|---|---|
| n >= 2 | Standard DiffBind path with biological replicates (publication-grade) |
| n = 1 | Auto-split BAM into 3 deterministic pseudo-replicates (`--n_pseudo_replicates`); fold-change estimates remain valid; p-values are optimistic -- **not for biological inference** |
| n = 0 (non-input) | Pipeline fails fast with a clear error |

## Parameters

Run `nextflow run <pipeline> --help` for the full list with descriptions.

Key parameters and their defaults:
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

## Citation

If you use this pipeline, please cite the preprint:
> Wang T., Jessa S., Marinov G.K., Klemm S., Kundaje A., Greenleaf W.J. (2025).
> Sensitive, direct detection of non-coding off-target base editor unwinding and editing in primary cells.
> bioRxiv. doi:10.1101/2025.09.25.678665

