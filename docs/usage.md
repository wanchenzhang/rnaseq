# computational-workflows/rnaseq: Usage

The complete instructions, parameter descriptions and troubleshooting notes are in the [README](../README.md).

## Run the integration test

Use Linux or WSL2 with Nextflow >=25.10.4 and a working Docker installation:

```bash
git clone https://github.com/wanchenzhang/rnaseq.git
cd rnaseq
nextflow run . -profile test,docker --outdir results_test
```

The profile supplies public FASTQ URLs and matching test references. It covers five samples with single-end, paired-end and multiple-run inputs. The supplied transcript FASTA skips GFFREAD. The small subsets are intended for integration testing. The `test_full` profile uses the same inputs and is not a full-size benchmark.

## Use your own inputs

Prepare a CSV samplesheet:

```csv
sample,fastq_1,fastq_2
SAMPLE_A,/absolute/path/SAMPLE_A_1.fastq.gz,/absolute/path/SAMPLE_A_2.fastq.gz
SAMPLE_B,/absolute/path/SAMPLE_B.fastq.gz,
```

Leave `fastq_2` empty for single-end reads. Repeat a sample ID for technical runs of the same sample; CAT_FASTQ merges them before trimming. Do not mix single-end and paired-end runs under one ID. Sample IDs cannot contain spaces. An [example layout](../assets/samplesheet.csv) is included; replace its placeholder paths with your own files.

```bash
nextflow run . \
    -profile docker \
    --input samplesheet.csv \
    --fasta /absolute/path/genome.fa \
    --gtf /absolute/path/genes.gtf \
    --outdir results
```

Genome FASTA and GTF must come from matching references. Add `--transcript_fasta /absolute/path/transcripts.fasta` to use matching precomputed transcripts; otherwise GFFREAD generates them. Review STAR index settings and resource requests in `conf/modules.config` before using another organism or larger datasets.

## Select preprocessing

Trim Galore is the default. Use `--trimmer fastp` to select fastp, or run:

```bash
nextflow run . \
    -profile test,docker,fastp_optional \
    --outdir results_test_fastp_optional
```

The optional profile enables polyX trimming and paired-end overlap correction. PolyG behavior is controlled with `--fastp_poly_g auto`, `on` or `off`. For custom Boolean settings, use unquoted `true` or `false` in a configuration profile. UMI-aware counting and overlapping-mate merging are not implemented.

## Check outputs

For the default Trim Galore test:

```bash
python bin/validate_outputs.py \
    --samplesheet https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/samplesheet/v3.10/samplesheet_test.csv \
    --outdir results_test \
    --trimmer trimgalore
```

For fastp, select its output directory and use `--trimmer fastp`. See [output documentation](output.md) and the README for the checks and their limitations.
