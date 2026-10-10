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

## Predicting gRNA Binding Sites

Set `--guidescan_grna` to enable these steps.

Given a base editor gRNA (protospacer + PAM), the pipeline reports the most likely
binding site per differential peak, via GuideScan2 and CFD scoring, in two stages:

**Stage 1 - `GUIDESCAN_SEED`** Greps each peak and its reverse
complement for an exact *seed* match (a substring of spacer + PAM: by default the
PAM-proximal 5 nt plus the PAM, so `GGCACTGCGGCTGGAGGTGG` + `NGG` -> `GGTGGNGG`).
Peaks with no hit are reported with `seed_hits = 0`.

**Stage 2 - `GUIDESCAN_ENUMERATE`** Runs GuideScan2 `enumerate` over the `seed_hits == 0` peaks. This step is computationally expensive for large differential peak sets: pass `--skip_guidescan_enumerate` to skip the second stage for a fast first pass of stage 1 alone.


```bash
nextflow run main.nf --stage diffbind -resume \
    --guidescan_grna CCCTACCTGTCACCAGGACC --guidescan_pam NGG \
    --samplesheet samples.csv --refdir /path/to/refs --outdir results/
```

A single `--guidescan_grna` applies globally to every case in the run. For multi-case runs that require separate binding sequences, write a `guidescan_params.json` keyed by `case_id`:
```
{
  "ABE8e_HEK4": { "grna": "GGCACTGCGGCTGGAGGTGG" },
  "ABE8e_CD7":  { "grna": "CCCTACCTGTCACCAGGACC" }
}
```
Other overridable parameters include: `pam`, `seed`, `mismatches`, `rna_bulges`,
`dna_bulges`, `chunk_size`, `summit_width`, `enumerate_mismatches`, `enumerate_chunk_size`

```bash
nextflow run main.nf --stage diffbind -resume \
    --guidescan_params guidescan_params.json \
    --samplesheet samples.csv --refdir /path/to/refs --outdir results/
```

## The Non-Differential Arm

`--nondiff_arm` reruns guidescan seed, variant calling and protospacer over
`nondiff_peaks.bed`, MACS2 treatment peaks that do not overlap a differential peak. Accounts
for highly transcribed loci with large endogenous ssDNA signal (un-caught gRNA-directed edits).

```bash
nextflow run main.nf --stage diffbind -resume --nondiff_arm \
    --guidescan_params guidescan_params.json \
    --samplesheet samples.csv --refdir /path/to/refs --outdir results/
```

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
    <case>/                                # differential_peaks.bed, nondiff_peaks.bed, full TSV, volcano, MA plot, dba_object.rds, method.txt
  07_variants/<case>/${sample_id}/         # *.vcf, *_no-alt-chr.vcf, *_CT.vcf, *_AG.vcf
  08_guidescan/<case>/                     # binding_sites_seed.tsv, binding_sites_enumerate.tsv
  09_variants_subtracted/<case>/           # <case>_vs_control_relcontrol*.vcf
  10_protospacer/                          # cohort figures + CSVs
  11_nondiff/                              # --nondiff_arm only
    <case>/guidescan/                      # binding_sites_seed.tsv, scored_peaks.bed
    <case>/variants/${sample_id}/          # *.vcf
    <case>/variants_subtracted/            # <case>_vs_control_relcontrol*.vcf
    protospacer/                           # cohort figures + CSVs, non-differential arm
  pipeline_info/                           # execution_report, timeline, trace, DAG SVG
```

## Pipeline Structure

The pipeline is split into two sub-workflows (selectable with `--stage`, see [Staged Execution](#staged-execution)):
- **`KAS_ANALYZER`** (stage 1) - trim/align through the per-case DiffBind samplesheets
- **`DIFFERENTIAL`** (stage 2) - differential peaks, variants, guidescan, protospacer

| Process | Stage | Source file | Description |
|---|---|---|---|
| `TRIM_ALIGN` | KAS_ANALYZER | `trim_align.nf` | trim_galore + KAS-Analyzer KAS-seq (Bowtie2, KAS-Analyzer's default aligner) |
| `PEAK_CALL` | KAS_ANALYZER | `analyze.nf` | MACS2 via KAS-Analyzer, treatment vs background |
| `KAS_STATISTICS` | KAS_ANALYZER | `analyze.nf` | Cohort RPM statistics -> ratios.txt |
| `NORMALIZE_BIGWIG` | KAS_ANALYZER | `analyze.nf` | RPM-ratio normalized bedGraph + UCSC bigWig |
| `QC_METRICS` | KAS_ANALYZER | `analyze.nf` | FRiP + multiBigwigSummary correlation |
| `DEDUP_CONSENSUS` | KAS_ANALYZER | `picard_gencore_dedup.nf` | picard MarkDuplicates + gencore consensus dedup |
| `SPLIT_BAM_PSEUDOREPS` | KAS_ANALYZER | `split_bam_pseudoreps.nf` | (conditional) deterministic BAM partition for n=1 conditions |
| `BUILD_DIFFBIND_SAMPLESHEET` | DIFFERENTIAL | `differential_peaks.nf` | Fan the cohort into one DiffBind samplesheet per treatment case (stage 1 stops here) |
| `DIFFERENTIAL_PEAKS` | DIFFERENTIAL | `differential_peaks.nf` | DiffBind (DESeq2 backend); one control-vs-case contrast per case, output under `06_differential/<case>/` |
| `CALL_VARIANTS` | DIFFERENTIAL | `mutect2_call_variants.nf` | Mutect2 tumor-only + germline + PoN, per case, restricted to that case's differential peaks |
| `GUIDESCAN_SEED` | DIFFERENTIAL | `guidescan.nf` | Exact-seed grep + bounded-DP CFD scoring; most-likely gRNA binding site per differential peak (requires `--guidescan_grna` or `--guidescan_params`) |
| `GUIDESCAN_ENUMERATE` | DIFFERENTIAL | `guidescan.nf` | GuideScan2 enumerate + CFD over remaining peaks (disable with `--skip_guidescan_enumerate`) |
| `VARIANT_SUBTRACT` | DIFFERENTIAL | `variant_subtract.nf` | Per case: merge control and treatment replicate VCFs, subtract control positions, AD/DP filter, split CT/AG |
| `PROTOSPACER` | DIFFERENTIAL | `protospacer.nf` | Cohort-level: variant position within the predicted protospacer, per substitution type (disable with `--skip_protospacer`) |
| `GUIDESCAN_SEED_NONDIFF` | DIFFERENTIAL | `guidescan.nf` | non-differential arm: seed grep over `nondiff_peaks.bed` (`--nondiff_arm`) |
| `CALL_VARIANTS_NONDIFF` | DIFFERENTIAL | `mutect2_call_variants.nf` | non-differential arm: Mutect2 over the arm's scored peaks |
| `VARIANT_SUBTRACT_NONDIFF` | DIFFERENTIAL | `variant_subtract.nf` | non-differential arm: treatment minus control |
| `PROTOSPACER_NONDIFF` | DIFFERENTIAL | `protospacer.nf` | non-differential arm: separate figure set, no on-target |

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
| `--fdr` | 0.05 | DiffBind FDR cutoff (cohort default) |
| `--lfc` | 1 | DiffBind log2 fold-change cutoff (cohort default) |
| `--min_overlap` | 2 | DiffBind consensus minOverlap (cohort default) |
| `--greylist` | false | Apply DiffBind greylist masking (cohort default) |
| `--n_pseudo_replicates` | 3 | Pseudo-rep count when n=1 |
| `--seed` | 42 | Global random seed |
| `--diffbind_params` | (none) | JSON file of per-case DiffBind overrides for stage 2 (see [Adjusting Differential Expression Parameters](#adjusting-differential-expression-parameters)) |
| `--diffbind_samplesheet_dir` | (none) | Directory of edited per-case samplesheets for stage 2 (see [Staged execution](#staged-execution)) |
| `--guidescan_grna` | (none) | gRNA spacer (20 nt, no PAM); guidescan runs only when this or `--guidescan_params` is set (see [Predicting gRNA Binding Sites](#predicting-grna-binding-sites)) |
| `--guidescan_params` | (none) | JSON file of per-case gRNA/seed overrides, for runs whose cases target different loci |
| `--guidescan_pam` | NGG | PAM pattern (N = wildcard) |
| `--guidescan_seed` | (derived) | Stage-1 seed; defaults to the PAM-proximal 5 nt of the spacer + the PAM |
| `--guidescan_mismatches` | 10 | Stage 1: max mismatches, applied as a post-scoring filter (not a search budget) |
| `--guidescan_rna_bulges` | 2 | Max 1-nt RNA bulges (both stages) |
| `--guidescan_dna_bulges` | 2 | Max 1-nt DNA bulges (both stages) |
| `--guidescan_chunk_size` | 150 | Stage 1: peaks per padded-extraction unit |
| `--guidescan_summit_width` | 0 | Stage 1: restrict the seed search to the central N bp of each peak; 0 = full peak (**leave at 0**, see the warning above) |
| `--skip_guidescan_enumerate` | false | Skip stage 2 (GuideScan2 enumerate), leaving stage 1's seed-anchored results only |
| `--guidescan_enumerate_mismatches` | 8 | Stage 2 search budget; separate from `--guidescan_mismatches` and deliberately lower |
| `--guidescan_enumerate_chunk_size` | 15 | Stage 2: peaks per GuideScan2 index+enumerate unit |
| `--skip_protospacer` | false | Skip the cohort-level protospacer position analysis (see [Protospacer Position Analysis](#protospacer-position-analysis)) |
| `--protospacer_palette` | (none) | JSON map of `case_id` -> hex colour, for exact figure colour parity |
| `--protospacer_both_strands` | false | `variant_altcount_in_protospacer` plots: also show TS edits below the axis (default: NTS only; the CSV always has both strands) |
| `--nondiff_arm` | false | Second pass over MACS2 peaks that were not called differential (see [The Non-Differential Arm](#the-non-differential-arm)) |

## Tool dependencies

Tool dependencies and pinned versions live in:
- [env/env.yaml](env/env.yaml) - conda environment spec (bioinformatics tools: BWA, samtools, bcftools, bedtools, MACS2, trim_galore, FastQC, deepTools, picard, gencore, GATK4, GuideScan2, UCSC bedGraphToBigWig).
- KAS-Analyzer 2.0.0 is cloned from `https://github.com/Ruitulyu/KAS-Analyzer.git`

## Citation

If you use this pipeline, please cite the preprint:
> Wang T., Jessa S., Marinov G.K., Klemm S., Kundaje A., Greenleaf W.J. (2025).
> Sensitive, direct detection of non-coding off-target base editor unwinding and editing in primary cells.
> bioRxiv. doi:10.1101/2025.09.25.678665
