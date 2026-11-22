# beCasKAS

Reproducible Nextflow pipeline implementing the beCasKAS computational analysis from
[Wang T. et al., bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.09.25.678665v1).
Detects non-coding off-target base editor unwinding and editing.

## Quick start
### Requirements
- [Nextflow](https://www.nextflow.io/) >= 25.04.0
- [Apptainer](https://apptainer.org/docs/admin/main/installation.html) >= 1.2.0
- ~50 GB scratch space + reference data (~30 GB for hg38 + GATK resources)

### Run Command

Clone the repository (check out a release tag), build the container (see [Container support](#container-support)), stage the references (see [Required References](#required-references)), then run:

```bash
git clone https://github.com/GreenleafLab/beCasKAS.git
cd beCasKAS
apptainer build becaskas.sif env/becaskas.def

nextflow run main.nf \
    --samplesheet samples.csv \
    --refdir /path/to/refs \
    --outdir results/ \
    --partition <your_partitions>
```

The pipeline does not download references per run. They must be staged in a directory and passed to `--refdir`.

Run `nextflow run main.nf --help` to list every parameter with its description and default.

By default, the pipeline is configured for a SLURM/HPC cluster (`-profile slurm`). Each task is submitted to `--partition`, whose default is `normal`, Sherlock's general partition: set it to partitions your account can use, especially on other clusters. To submit a thin driver job see [submit_example.sbatch](examples/submit_example.sbatch) for a template. To run on a single local machine, add `-profile local`, which switches to the local executor and automatically caps each job's resource request to the CPUs/RAM available on that machine:

```bash
nextflow run main.nf -profile local \
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

Apptainer is the default container runtime (configured in [nextflow.config](nextflow.config)) and works for both the default SLURM/HPC target and `-profile local`. A Docker profile is defined, but no Docker image is published. Build your own following [env/becaskas.def](env/becaskas.def) and pass it with `--container`.

## Required References

Supply reference genome for read mapping, GATK variant calling resource files. Stage the following reference layout under `--refdir`:

```
hg38.fasta                              # UCSC hg38 fasta
hg38.fasta.fai                          # samtools faidx
hg38.dict                               # samtools dict

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

All reference files live flat in `${refdir}/`. The pipeline expects exactly these filenames. `--genome` sets the prefix: `--genome hg38_no_alt` reads `hg38_no_alt.*` instead of `hg38.*`.

The staging steps below use the `becaskas.sif` container, which already bundles `samtools` and `bowtie2-build`. Run them from inside `${refdir}`.

1. Download UCSC `hg38.fa.gz` from [UCSC HG38 Download Page](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz), `gunzip` and rename to `hg38.fasta`.

2. Build the fasta index/dict and the Bowtie2 index (~30 to 60 minutes; produces the six `hg38.*.bt2` files). Point to your `becaskas.sif`:
```bash
apptainer exec becaskas.sif samtools faidx hg38.fasta
apptainer exec becaskas.sif samtools dict hg38.fasta -o hg38.dict
apptainer exec becaskas.sif bowtie2-build hg38.fasta hg38
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

Declares sample libraries for analysis (editor + gRNA, editor, background)

CSV with columns: `sample_id`, `fastq_1`, `fastq_2`, `condition`, `replicate`, `background_id`. See [samplesheet_example.csv](examples/samplesheet_example.csv).

| Column | Required | Description |
|---|---|---|
| `sample_id` | yes | Unique sample identifier (letters, digits, `_`, `.`, `-`) |
| `fastq_1` | yes | Path/URL to R1 FASTQ (`.fastq`, `.fq`, optionally `.gz`) |
| `fastq_2` | yes | Path/URL to R2 FASTQ. The pipeline is paired-end only. |
| `condition` | yes | Letters, digits and `_`. Reserved: `background` (MACS2 input) and `control` (shared reference). Any other label defines a treatment case compared against `control`; the label becomes the `case_id` used in output paths and JSON overrides. |
| `replicate` | yes | Replicate number within condition (positive integer) |
| `background_id` | conditional | Required on every non-`background` row: the sample_id of the `background` row used as MACS2 `-c` |

## Staged Execution


Splits run between per-sample ssDNA measurements and cross-sample comparisons. The pipeline runs as two stages: `KAS_ANALYZER` (trim/align through the per-case DiffBind samplesheets) and `DIFFERENTIAL` (differential peaks, variants, guidescan, protospacer, plus the optional non-differential arms).

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

Determines and separates out gRNA-directed unwinding sites (differential peaks) compared to editor-only controls. 

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

Aligns the gRNA within a differential (unwound) peak, locates the editor bound protospacer and measures match to the on-target sequence.

Set `--guidescan_grna` (or `--guidescan_params`) to enable these steps.

Given a base editor gRNA (protospacer + PAM), the pipeline reports the most likely
binding site per differential peak, scored with Doench 2016 CFD, in two stages:

**Stage 1 - `GUIDESCAN_SEED`** Greps each peak and its reverse
complement for an exact *seed* match (a substring of spacer + PAM: by default the
PAM-proximal 5 nt plus the PAM, so `CCCTACCTGTCACCAGGACC` + `NGG` -> `GGACCNGG`).
Peaks with no hit are reported with `seed_hits = 0`.

**Stage 2 - `GUIDESCAN_ENUMERATE`** Runs GuideScan2 `enumerate` over the `seed_hits == 0` peaks. This step is computationally expensive for large differential peak sets: pass `--skip_guidescan_enumerate` to skip the second stage for a fast first pass of stage 1 alone.


```bash
nextflow run main.nf --stage diffbind -resume \
    --guidescan_grna CCCTACCTGTCACCAGGACC --guidescan_pam NGG \
    --samplesheet samples.csv --refdir /path/to/refs --outdir results/
```

A single `--guidescan_grna` applies globally to every case in the run. To set the gRNA per case, write a `guidescan_params.json` keyed by `case_id`.
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

## The Non-Differential Arms

Looks for editing at gRNA-bound sites in high endogenous ssDNA regions that may not be caught by the differential arm due to high existing signal.

Non-differential arms are off by default and run independent workflows to each other. They are enabled with the flags outlined below and both must be specified to run them together.

### Sharp / Narrow (`--sharp_nondiff`)

Covers gRNA-bound sites that sit in concentrated heavily transcribed loci (high endogenous ssDNA signal).

`--sharp_nondiff` runs guidescan seed, variant calling and protospacer over
`nondiff_peaks.bed`, MACS2 treatment peaks that do not overlap a differential peak. Accounts for highly transcribed loci with large endogenous ssDNA signal (un-caught gRNA-directed edits).

```bash
nextflow run main.nf --stage diffbind -resume --sharp_nondiff \
    --guidescan_params guidescan_params.json \
    --samplesheet samples.csv --refdir /path/to/refs --outdir results/
```

### Broad (`--broad_nondiff`)

Covers gRNA-bound sites inside broad ssDNA domains, where the unwinding is spread across a long stretch of chromatin.

`--broad_nondiff` runs EPIC2 to call broad peaks, then guidescan seed, variant calling and protospacer over `broad_regions.bed`. Finds gRNA-directed off-targets in broad KAS-seq domains.

We recommend `--genome hg38_no_alt` as EPIC2 broad peaks tend to map to alt contigs at a higher rate.

**NOTE:** Broad regions overlapping sharp peaks are not subtracted by default. A site inside both may be scored by both arms, but `12_broad_nondiff/<case>/arm_overlap.tsv` and `arm_overlap_venn.pdf` will report the overlaps.

```bash
nextflow run main.nf --stage diffbind -resume --broad_nondiff \
    --guidescan_params guidescan_params.json \
    --samplesheet samples.csv --refdir /path/to/refs --outdir results/
```

## Output

Outputs are separated and numbered by pipeline sub-process outputs:

```
results/
  01_trim_align/${sample_id}/              # sorted BAM + index, .bed, .bg, mapping summary, FastQC + trimming reports
  02_peaks/${sample_id}/                   # MACS2 peaks (.bed, .xls, .narrowPeak, summits) (+ *_input_epic2_peaks.bed with --broad_nondiff)
  03_normalization/                        # KAS_statistics_summary*.txt, ratios.txt (RPM scaling factors)
  03_bigwigs/${sample_id}/                 # *.nor.bigWig, *.nor.bg
  04_qc/                                   # FRiP report, bigWig correlation (all-deeptools.pdf, all.tab, all.npz)
  05_dedup/${sample_id}/                   # *_picard-gencore-dedup.bam + index, picard metrics, gencore report
    pseudoreps/                            # (only if n=1) pseudo-rep BAMs
  06_differential/
    samplesheets/                          # diffbind_samplesheet_<case>.csv (editable; --diffbind_samplesheet_dir)
    <case>/                                # differential_peaks.bed, nondiff_peaks.bed, differential_peaks_full.tsv, volcano.pdf,
                                           #   ma_plot.pdf, dba_object.rds, differential_peaks_method.txt
  07_variants/<case>/${sample_id}/         # Mutect2 VCFs: *.vcf (raw), *_no-alt-chr.vcf, *_CT.vcf, *_AG.vcf
  08_guidescan/<case>/                     # binding_sites_seed.tsv, binding_sites_enumerate.tsv, binding_sites_seed_raw_hits.tsv,
                                           #   scored_peaks.bed, no_seed_peaks.bed (stage 2 input)
  09_variants_subtracted/<case>/           # <case>_vs_control_relcontrol_filt_{no-alt-chr,CT,AG}.vcf
  10_protospacer/                          # cohort figures + CSVs (edit position within the protospacer), vcf/
  11_sharp_nondiff/                        # --sharp_nondiff only
    <case>/guidescan/                      # binding_sites_seed.tsv, scored_peaks.bed
    <case>/variants/${sample_id}/          # *.vcf
    <case>/variants_subtracted/            # <case>_vs_control_relcontrol_filt_{no-alt-chr,CT,AG}.vcf
    protospacer/                           # cohort figures + CSVs, non-differential arm
  12_broad_nondiff/                        # --broad_nondiff only
    <case>/                                # broad_regions.bed, broad_regions_summary.txt, arm_overlap.tsv, arm_overlap_venn.pdf
    <case>/guidescan/                      # binding_sites_seed.tsv (all anchors), scored_peaks.bed
    <case>/variants/${sample_id}/          # *.vcf
    <case>/variants_subtracted/            # <case>_vs_control_relcontrol_filt_{no-alt-chr,CT,AG}.vcf
    protospacer/                           # cohort figures + CSVs, broad arm
  pipeline_info/                           # report.html, timeline.html, trace.txt, pipeline_dag.svg (.dot without Graphviz)
  pipeline.success                         # or pipeline.failure when the run ends (pipeline.running while it runs)
```

## DiffBind Replicate Handling

Allows ssDNA comparison when a condition has only one library. Splits library's reads into pseudo-replicates to stand in for replication.

| Samples per condition | Behavior |
|---|---|
| n >= 2 | Standard DiffBind path with biological replicates |
| n = 1 | Auto-split BAM into 3 deterministic pseudo-replicates (`--n_pseudo_replicates`) |
| no `control` or no treatment rows | `BUILD_DIFFBIND_SAMPLESHEET` fails with an error naming the condition labels it found |

A case with pseudo-replicates on either side uses library-size normalization (`DBA_NORM_LIB`) and a stricter default log2FC cutoff of 2. `--lfc`.

## Parameters

Run `nextflow run main.nf --help` for the full list with descriptions.

Key parameters and their defaults:
| Param | Default | Role |
|---|---|---|
| `--samplesheet` | (required) | Samplesheet CSV (see [Samplesheet](#samplesheet)) |
| `--refdir` | (required) | Directory of staged references (see [Required References](#required-references)) |
| `--outdir` | results | Output directory |
| `--container` | `<repo>/becaskas.sif` | Container image used by every process (see [Container support](#container-support)) |
| `--stage` | all | Stage to run: `all` / `kas_analyzer` / `diffbind` |
| `--genome` | hg38 | Reference file prefix in `--refdir`, e.g. `hg38_no_alt` (see [Required References](#required-references)) |
| `--clip_r2` | 12 | Bases clipped from 5' of R2 (adaptase library prep) |
| `--fdr` | 0.05 | DiffBind FDR cutoff (cohort default) |
| `--lfc` | 1 (2 for pseudo-reps) | DiffBind log2 fold-change cutoff (cohort default) |
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
| `--guidescan_summit_width` | 0 | Stage 1: restrict the seed search to the central N bp of each peak; 0 = full peak |
| `--skip_guidescan_enumerate` | false | Skip stage 2 (GuideScan2 enumerate), leaving stage 1's seed-anchored results only |
| `--guidescan_enumerate_mismatches` | 8 | Stage 2 search budget; separate from `--guidescan_mismatches` and deliberately lower |
| `--guidescan_enumerate_chunk_size` | 15 | Stage 2: peaks per GuideScan2 index+enumerate unit |
| `--skip_protospacer` | false | Skip the cohort-level protospacer position analysis |
| `--protospacer_palette` | (none) | JSON map of `case_id` -> hex colour, for exact figure colour parity |
| `--protospacer_both_strands` | false | `variant_altcount_in_protospacer` plots: also show TS edits below the axis (default: NTS only; the CSV always has both strands) |
| `--sharp_nondiff` | false | Run arm for non-differential MACS2 sharp peaks (see [The Non-Differential Arms](#the-non-differential-arms)) |
| `--broad_nondiff` | false | Run arm for non-differential EPIC2 broad peaks (see [The Non-Differential Arms](#the-non-differential-arms))|
| `--partition` | normal | SLURM partition(s) for task jobs; set to partitions your account can use |
| `--clusterOptions` | (none) | Extra flags appended to every task's `sbatch` call, e.g. `--qos=long` |
| `--max_cpus` / `--max_memory` / `--max_time` | 16 / 128.GB / 48.h | Upper limits on any one task's request |

## Tool dependencies

Tool dependencies and pinned versions live in:
- [env/env.yaml](env/env.yaml) - conda environment spec: trim_galore, cutadapt, FastQC, Bowtie2, samtools, bcftools, bedtools, MACS2, EPIC2, deepTools, picard, gencore, GATK4, GuideScan2, UCSC bedGraphToBigWig, and R with DiffBind, DESeq2 and edgeR.
- [env/becaskas.def](env/becaskas.def) - Apptainer recipe; also clones KAS-Analyzer v2.0.0 from `https://github.com/Ruitulyu/KAS-Analyzer.git`

## Citation

If you use this pipeline, please cite the preprint:
> Wang T., Jessa S., Liu E.S.X., Marinov G.K., Klemm S., Russkamp N.F., Porteus M.H., Kundaje A., Greenleaf W.J. (2025).
> Sensitive, direct detection of non-coding off-target base editor unwinding and editing in primary cells.
> bioRxiv. doi:10.1101/2025.09.25.678665
