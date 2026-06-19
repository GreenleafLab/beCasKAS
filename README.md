# beCasKAS

Reproducible Nextflow pipeline implementing the beCasKAS computational analysis from
[Wang T. et al., bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.09.25.678665v1).
Detects non-coding off-target base editor unwinding and editing in primary cells.

## Quick start

### Requirements
- [Nextflow](https://www.nextflow.io/) >= 25.10.0
- [Apptainer](https://apptainer.org/docs/admin/main/installation.html) >= 1.2.0
- ~50 GB scratch space + reference data (~30 GB for hg38 + GATK resources)

### Run Command

```bash
nextflow run GreenleafLab/beCasKAS \
    --samplesheet samples.csv \
    --refdir /path/to/refs \
    --outdir results/
```

The pipeline does not download references per run. They must be staged in a directory and passed to `--refdir`.

By default, the pipeline is configured for a SLURM/HPC cluster (`-profile slurm`). To submit a thin driver job see [submit_example.sbatch](examples/submit_example.sbatch) for a template. To run on a single local machine, add `-profile local`, which switches to the local executor and automatically caps each job's resource request to the CPUs/RAM available on that machine:

```bash
nextflow run GreenleafLab/beCasKAS -profile local \
    --samplesheet samples.csv \
    --refdir /path/to/refs \
    --outdir results/
```

## Container support

The pipeline ships an Apptainer recipe under [env/becaskas.def](env/becaskas.def). Build the SIF locally from the repository root:

```bash
apptainer build --fakeroot becaskas.sif env/becaskas.def
```

- `--fakeroot` enables unprivileged builds and works on most modern Linux distributions (requires user namespace support). HPCs that disable user namespaces will require the admin to build the SIF and distribute it.
- HPCs may have login/interactive download and memory limits. It is typically advised to run a `sbatch` SLURM job to build the container (See [container_build_example.sbatch](examples/container_build_example.sbatch))

By default the pipeline expects `becaskas.sif` at the repository root. To use a SIF in a different location, override the default:

```bash
nextflow run main.nf --container /abs/path/to/becaskas.sif ...
```

Apptainer is the default container runtime (configured in [nextflow.config](nextflow.config)) and works for both the default SLURM/HPC target and `-profile local`. A Docker profile is pre-configured, but requires separate installation.

## Required References

Stage the following reference layout under `--refdir`:

```
hg38.fasta                              # UCSC hg38 fasta
hg38.fasta.fai                          # samtools faidx
hg38.dict                               # samtools dict
hg38.genome                             # chromosome sizes (UCSC chrom.sizes)

hg38.1.bt2                              # bowtie2 indices
hg38.2.bt2
hg38.3.bt2
hg38.4.bt2
hg38.rev.1.bt2
hg38.rev.2.bt2

1000g_pon.hg38.vcf.gz                   # GATK resources
1000g_pon.hg38.vcf.gz.tbi
af-only-gnomad.hg38.vcf.gz
af-only-gnomad.hg38.vcf.gz.tbi
```

All reference files live flat in `${refdir}/`. The pipeline expects exactly these filenames.

Suggesting staging steps directly use the generated `becaskas.sif` container as dependencies for `samtools` and `bowtie2-build` commands are already bundled inside. Run the following from inside `${refdir}`.

1. Download UCSC `hg38.fa.gz` from [UCSC HG38 Download Page](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz),`gunzip` and rename to `hg38.fasta`.

2. Build the fasta index/dict, the Bowtie2 index (~30 to 60 minutes; produces the six `hg38.*.bt2` files), and the chrom.sizes (specify your specific path to `becaskas.sif`):
```bash
apptainer exec becaskas.sif samtools faidx hg38.fasta
apptainer exec becaskas.sif samtools dict hg38.fasta -o hg38.dict
apptainer exec becaskas.sif bowtie2-build hg38.fasta hg38
cut -f1,2 hg38.fasta.fai > hg38.genome   # chrom.sizes from the .fai above, or download from UCSC
```

3. Download the GATK somatic resources:
```bash
BASE=https://storage.googleapis.com/gatk-best-practices/somatic-hg38
wget $BASE/af-only-gnomad.hg38.vcf.gz
wget $BASE/af-only-gnomad.hg38.vcf.gz.tbi
wget $BASE/1000g_pon.hg38.vcf.gz
wget $BASE/1000g_pon.hg38.vcf.gz.tbi
```

## Samplesheet

CSV with columns: `sample_id`, `fastq_1`, `fastq_2`, `condition`, `replicate`, `background_id`. See [samplesheet_example.csv](examples/samplesheet_example.csv).

| Column | Required | Description |
|---|---|---|
| `sample_id` | yes | Unique sample identifier |
| `fastq_1` | yes | Path/URL to R1 FASTQ |
| `fastq_2` | no | Path/URL to R2 FASTQ (omit for single-end) |
| `condition` | yes | Reserved: `background` (MACS2 input) and `control` (shared reference). Any other label defines a treatment case compared against `control`. |
| `replicate` | yes | Replicate number within condition |
| `background_id` | conditional | sample_id of the `background` row used as MACS2 `-c` |

## Staged Execution

The pipeline runs in two stages, selected with `--stage`:

- `--stage all` (default) - full pipeline end to end.
- `--stage kas_analyzer` - run everything through the per-case DiffBind samplesheet generation, then stop. The samplesheets are published to `results/06_differential/samplesheets/` and the case summary is printed to the run log.
- `--stage diffbind` - run only the differential + variant stage. Use it with `-resume` after a `kas_analyzer` run so the trim/align/dedup/peak steps are served from cache:
```
nextflow run main.nf --stage kas_analyzer <args>

# inspect results/06_differential/samplesheets/ and the printed case summary
nextflow run main.nf -resume --stage diffbind <args>
```
`-resume` is required for `--stage diffbind` to reuse the cached upstream; without it the upstream re-runs from scratch. Because this is a resumed continuation, the earlier stage `work/` directories must still exist.

To edit which samples or paths go into a case, edit the published `diffbind_samplesheet_<case>.csv` files and point the second stage at its directory. (**NOTE**: modifying the `work/` or `results/` copy without passing the argument will not apply your edits, it must be passed via `--diffbind_samplesheet_dir`):
```
nextflow run main.nf -resume --stage diffbind --diffbind_samplesheet_dir results/06_differential/samplesheets <args>
```

### Adjusting Differential Expression Parameters

By default, multiple treatment cases (distinct `condition` labels) are each compared against the single shared `control`. Per-case differential parameters can be modified in the second stage (`--stage diffbind`) without re-running the rest of the pipeline. Write a `diffbind_params.json` file that maps each case to the parameters you would like to change:
```
    {
      "<case_id>": {
        "fdr": 0.01, "lfc": 2,
        "dba.count":     { "minOverlap": 1,
                           "summits": 100,
                           "score": "DBA_SCORE_READS" },
        "dba.analyze":   { "method": "DBA_EDGER" },
        "dba.normalize": { "normalize": "DBA_NORM_LIB" }
      }
    }
```
then run the differential stage, passing the file with `--diffbind_params`:
```
    nextflow run main.nf --stage diffbind -resume --diffbind_params diffbind_params.json ...
```
Only the overridden case's `DIFFERENTIAL_PEAKS` and its downstream will re-execute. All other cases and every upstream stage will resume from cache.

## Output

```
results/
  01_trim_align/${sample_id}/             # sorted.bam, .bed, .bg, fastqc reports
  02_peaks/${sample_id}/                   # MACS2 peak BED + xls
  03_bigwigs/${sample_id}/                 # *.nor.bigWig, *.nor.bg
  04_qc/                                   # FRiP report, correlation PDF/tab
  05_dedup/${sample_id}/                   # *_picard-gencore-dedup.bam
    pseudoreps/                            # (only if n=1) pseudo-rep BAMs
  06_differential/
    samplesheets/                          # diffbind_samplesheet_<case>.csv (editable; --diffbind_samplesheet_dir)
    <case>/                                # differential_peaks.bed, full TSV, volcano, dba_object.rds, method.txt
  07_variants/<case>/${sample_id}/         # *.vcf, *_no-alt-chr.vcf, *_CT.vcf, *_AG.vcf
  pipeline_info/                           # execution_report, timeline, trace, DAG SVG
```

## Pipeline Structure

The pipeline is split into two sub-workflows (selectable with `--stage`, see [Staged Execution](#staged-execution)):
- **`KAS_ANALYZER`** (stage 1) - steps 1-8
- **`DIFFERENTIAL`** (stage 2) - steps 9-10

| # | Process | Workflow | Source file | Description |
|---|---|---|---|---|
| 1 | `TRIM_ALIGN` | KAS_ANALYZER | `trim_align.nf` | trim_galore + KAS-Analyzer KAS-seq (Bowtie2, KAS-Analyzer's default aligner) |
| 2 | `PEAK_CALL` | KAS_ANALYZER | `analyze.nf` | MACS2 via KAS-Analyzer (`--macs2_mode`), treatment vs background |
| 3 | `KAS_STATISTICS` | KAS_ANALYZER | `analyze.nf` | Cohort RPM statistics -> ratios.txt |
| 4 | `NORMALIZE_BIGWIG` | KAS_ANALYZER | `analyze.nf` | RPM-ratio normalized bedGraph + UCSC bigWig |
| 5 | `QC_METRICS` | KAS_ANALYZER | `analyze.nf` | FRiP + multiBigwigSummary correlation |
| 6 | `DEDUP_CONSENSUS` | KAS_ANALYZER | `picard_gencore_dedup.nf` | picard MarkDuplicates + gencore consensus dedup |
| 7 | `SPLIT_BAM_PSEUDOREPS` | KAS_ANALYZER | `split_bam_pseudoreps.nf` | (conditional) deterministic BAM partition for n=1 conditions |
| 8 | `BUILD_DIFFBIND_SAMPLESHEET` | KAS_ANALYZER | `differential_peaks.nf` | Fan the cohort into one DiffBind samplesheet per treatment case (stage 1 stops here) |
| 9 | `DIFFERENTIAL_PEAKS` | DIFFERENTIAL | `differential_peaks.nf` | DiffBind (DESeq2 backend); one control-vs-case contrast per case, output under `06_differential/<case>/` |
| 10 | `CALL_VARIANTS` | DIFFERENTIAL | `mutect2_call_variants.nf` | Mutect2 tumor-only + germline + PoN, per case, restricted to that case's differential peaks |

## Replicate Handling

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
| `--stage` | all | Stage to run: `all` / `kas_analyzer` / `diffbind` |
| `--clip_r2` | 12 | Bases clipped from 5' of R2 (adaptase library prep) |
| `--macs2_mode` | sharp | MACS2 peak-calling mode |
| `--fdr` | 0.05 | DiffBind FDR cutoff (cohort default) |
| `--lfc` | 1 | DiffBind log2 fold-change cutoff (cohort default) |
| `--min_overlap` | 2 | DiffBind consensus minOverlap (cohort default) |
| `--greylist` | false | Apply DiffBind greylist masking (cohort default) |
| `--n_pseudo_replicates` | 3 | Pseudo-rep count when n=1 |
| `--seed` | 42 | Global random seed |
| `--diffbind_params` | (none) | JSON file of per-case DiffBind overrides for stage 2 (see [Adjusting Differential Expression Parameters](#adjusting-differential-expression-parameters)) |
| `--diffbind_samplesheet_dir` | (none) | Directory of edited per-case samplesheets for stage 2 (see [Staged execution](#staged-execution)) |

## Tool dependencies

Tool dependencies and pinned versions live in:
- [env/env.yaml](env/env.yaml) - conda environment spec (bioinformatics tools: BWA, samtools, bcftools, bedtools, MACS2, trim_galore, FastQC, deepTools, picard, gencore, GATK4, UCSC bedGraphToBigWig).
- KAS-Analyzer 2.0.0 is cloned from `https://github.com/Ruitulyu/KAS-Analyzer.git`

## Citation

If you use this pipeline, please cite the preprint:
> Wang T., Jessa S., Marinov G.K., Klemm S., Kundaje A., Greenleaf W.J. (2025).
> Sensitive, direct detection of non-coding off-target base editor unwinding and editing in primary cells.
> bioRxiv. doi:10.1101/2025.09.25.678665
