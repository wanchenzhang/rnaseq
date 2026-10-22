# Development notes

This is a two-person course project. The main branch is `master`.

When changing the pipeline:

1. Connect modules in `workflows/rnaseq.nf` and set tool options in `conf/modules.config`.
2. Add new parameters to `nextflow.config` and `nextflow_schema.json`.
3. Update the README, output notes and tool references when needed.
4. Run the test and check its outputs:

```bash
nextflow run . -profile test,docker --outdir results_test
python bin/validate_outputs.py \
    --samplesheet https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/samplesheet/v3.10/samplesheet_test.csv \
    --outdir results_test \
    --trimmer trimgalore
```

If trimming changes, test the fastp branch too:

```bash
nextflow run . -profile test,docker,fastp_optional --outdir results_test_fastp_optional
```

Run code and formatting checks before committing:

```bash
prek run --all-files
nf-core pipelines lint
```

Some hooks fix files. If that happens, review the changes and run the checks again. GitHub runs lint checks after a push to `master`. A green lint result does not mean the full pipeline was run in CI.

Keep downloaded data, results, `work/`, caches and credentials out of Git. Keep the original license and credit when reusing nf-core modules. Do not move an existing release tag to a different commit.
