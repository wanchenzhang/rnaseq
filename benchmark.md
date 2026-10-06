# Trimming benchmark: fastp vs Trim Galore

This preliminary benchmark compares the two trimming branches of our custom RNA-seq pipeline using the same three-sample yeast test dataset.

## Dataset and workflow

- Organism: Saccharomyces cerevisiae.
- Samples: WT_REP1, WT_REP2 and RAP1_IAA_30M_REP1.
- Input: 50,000 paired-end read pairs per sample; 150,000 pairs in total.
- Raw read length: 101 bp.
- Both trimming branches feed trimmed paired-end reads into the same downstream FastQC, STAR and Salmon workflow.
- These small datasets are used for pipeline testing, rather than biological inference.

## Results

| Metric | fastp | Trim Galore |
| --- | ---: | ---: |
| Overall read-pair retention | 95.60% | 99.53% |
| Mean trimming task execution time | 1.967 s | 0.391 s |
| Peak memory range recorded in the Nextflow trace | 353.3 MB - 1.1 GB | 8.7 - 8.8 MB |
| STAR accepted alignment rate | 55.41% - 60.58% | 55.32% - 60.42% |
| Post-trimming per-base sequence quality | PASS in all six reports | PASS in all six reports |
| Post-trimming adapter-content check | PASS in all six reports | PASS in all six reports |
| Post-trimming overrepresented-sequence check | WARN in all six reports | WARN in all six reports |

Each branch produces six post-trimming FastQC reports: three samples multiplied by R1 and R2. These are not six biological samples.

## Interpretation

- Trim Galore retained more read pairs: the retention difference was 3.93 percentage points.
- Trim Galore had a shorter mean trimming task execution time in these measurements. This does not establish that it is faster for larger datasets or other environments.
- The trace recorded lower peak memory values for Trim Galore. These are observed measurements, not the memory requested in the process configuration.
- The STAR accepted alignment-rate ranges were similar. Range summaries alone do not establish a meaningful improvement for individual samples.
- Post-trimming FastQC statuses were identical for the quality, adapter-content and overrepresented-sequence checks. These categorical statuses do not identify a better trimming tool or demonstrate that the underlying quantitative metrics are identical.

## Measurement scope and limitations

The values above reproduce the preliminary measurements collected for this project. They have not been independently recalculated in this document.

Read-pair retention is the total number of retained pairs divided by the total input pairs. STAR accepted alignment includes unique alignments and permitted multiple alignments, with trimmed input pairs as its denominator. It should not be directly compared with Salmon's alignment-mode mapped percentage, which uses a different input denominator.

Timing refers to trimming task execution, not the complete pipeline's wall-clock runtime. Comparisons should use the same trace field, hardware, software versions and resource settings, and should exclude cached tasks. Short tasks are sensitive to overhead and memory-sampling resolution; the unusually low Trim Galore memory values should be checked against the original trace field and process measurement before making broader performance claims.

Tool defaults and filtering rules can differ. Lower retention alone does not demonstrate worse trimming, and a FastQC PASS does not mean that adapter content is exactly zero. Biological conclusions should not be drawn from these small test subsets.

For reproducibility, retain the original execution traces, commands, module configuration, tool versions and MultiQC reports alongside this comparison.
