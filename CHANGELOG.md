# computational-workflows/rnaseq: Changelog

## Unreleased

- Exclude the report PDF from the template-placeholder check while retaining checks on source files.
- Correct repository metadata and update usage and output documentation.
- Replace the unused viral-reconstruction full-test example with an explicitly documented alias of the working integration test.
- Normalize formatting and ignore local nf-test caches.

## v1.0.0 - 2026-10-09

- Develop a short-read bulk RNA-seq workflow from the nf-core template and reusable modules.
- Support single-end, paired-end and multiple sequencing runs per sample through CAT_FASTQ.
- Provide Trim Galore and fastp branches, with optional fastp preprocessing features.
- Connect STAR alignment and internal sorting, Picard duplicate marking, Salmon quantification and tximport aggregation.
- Generate MultiQC reports and Nextflow execution records.
- Add a five-sample public test profile and Python output validation.
