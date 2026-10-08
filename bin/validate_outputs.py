#!/usr/bin/env python3
"""Validate single/paired-end, multi-run RNA-seq outputs using the standard library."""
import argparse
import csv
import gzip
import math
import io
from urllib.request import urlopen
from itertools import zip_longest
from pathlib import Path


def load_samples(path):
    """Group runs by sample ID; FASTQ input files are not downloaded."""
    if str(path).startswith(('https://', 'http://')):
        with urlopen(str(path), timeout=30) as response:
            text = response.read().decode('utf-8-sig')
    else:
        text = Path(path).read_text(encoding='utf-8-sig')
    reader = csv.DictReader(io.StringIO(text, newline=''))
    if not {'sample', 'fastq_1', 'fastq_2'} <= set(reader.fieldnames or []):
        raise ValueError('Samplesheet requires sample, fastq_1 and fastq_2 columns')
    samples = {}
    seen = set()
    for line, row in enumerate(reader, 2):
        if None in row or any(row.get(k) is None for k in ('sample', 'fastq_1', 'fastq_2')):
            raise ValueError(f'Samplesheet row {line}: incorrect column count')
        name = row['sample'].strip()
        if not name or any(c.isspace() for c in name) or '/' in name or '\\' in name or name in ('.', '..'):
            raise ValueError('Sample IDs must be nonempty, without spaces or path separators')
        r1, r2 = row['fastq_1'].strip(), row['fastq_2'].strip()
        if not r1:
            raise ValueError(f'{name}: fastq_1 is required')
        run = (name, r1, r2)
        if run in seen:
            raise ValueError(f'{name}: duplicate sequencing-run row')
        seen.add(run)
        paired = bool(r2)
        if name in samples and samples[name]['paired'] != paired:
            raise ValueError(f'{name}: cannot mix single-end and paired-end runs')
        entry = samples.setdefault(name, {'paired': paired, 'runs': 0})
        entry['runs'] += 1
    if not samples:
        raise ValueError('Samplesheet is empty')
    return samples


def fastq_ids(path, mate=None):
    with gzip.open(path, 'rt', encoding='ascii') as handle:
        record = 0
        while True:
            header = handle.readline()
            if not header:
                return
            record += 1
            seq = handle.readline().rstrip('\r\n')
            plus = handle.readline().rstrip('\r\n')
            quality = handle.readline().rstrip('\r\n')
            if not header.startswith('@') or not plus.startswith('+'):
                raise ValueError(f'{path.name}: invalid four-line FASTQ record {record}')
            if not seq or len(seq) != len(quality):
                raise ValueError(f'{path.name}: sequence/quality length mismatch at record {record}')
            if any(ord(c) < 33 or ord(c) > 126 for c in quality):
                raise ValueError(f'{path.name}: invalid quality character at record {record}')
            fields = header[1:].strip().split()
            if not fields:
                raise ValueError(f'{path.name}: empty read ID at record {record}')
            read_id = fields[0]
            if read_id.endswith(('/1', '/2')):
                if mate is not None and not read_id.endswith(f'/{mate}'):
                    raise ValueError(f'{path.name}: incorrect mate suffix at record {record}')
                read_id = read_id[:-2]
            if mate is not None and len(fields) > 1 and fields[1].startswith(('1:', '2:')):
                if not fields[1].startswith(f'{mate}:'):
                    raise ValueError(f'{path.name}: incorrect mate field at record {record}')
            yield read_id


def check_pair(r1, r2):
    count = 0
    for count, (id1, id2) in enumerate(zip_longest(fastq_ids(r1, 1), fastq_ids(r2, 2)), 1):
        if id1 is None or id2 is None:
            raise ValueError(f'R1/R2 record counts differ at pair {count}')
        if id1 != id2:
            raise ValueError(f'R1/R2 read IDs differ at pair {count}: {id1} / {id2}')
    if count == 0:
        raise ValueError('No read pairs remain')
    return count


def check_single(path):
    count = sum(1 for _ in fastq_ids(path))
    if not count:
        raise ValueError('No reads remain')
    return count


def output_fastqs(outdir, sample, paired, stage):
    if stage == 'cat':
        filenames = ([f'{sample}_1.merged.fastq.gz', f'{sample}_2.merged.fastq.gz']
                     if paired else [f'{sample}.merged.fastq.gz'])
    elif stage == 'fastp':
        filenames = ([f'{sample}_R1.fastp.fastq.gz', f'{sample}_R2.fastp.fastq.gz']
                     if paired else [f'{sample}.fastp.fastq.gz'])
    else:
        filenames = ([f'{sample}_1_val_1.fq.gz', f'{sample}_2_val_2.fq.gz']
                     if paired else [f'{sample}_trimmed.fq.gz'])
    return [outdir / stage / filename for filename in filenames]


def check_tpm(path, samples):
    with path.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader, [])
        if len(header) != len(set(header)) or not header or header[0] != 'gene_id':
            raise ValueError('Invalid or duplicate TPM column names; first column must be gene_id')
        observed = [name for name in header if name not in ('gene_id', 'gene_name')]
        if set(observed) != set(samples):
            raise ValueError(f'TPM sample columns do not match samplesheet: {observed}')
        indexes = {name: header.index(name) for name in samples}
        totals = {name: 0.0 for name in samples}
        genes = set()
        for line, row in enumerate(reader, 2):
            if len(row) != len(header):
                raise ValueError(f'TPM row {line}: incorrect column count')
            gene = row[0].strip()
            if not gene or gene in genes:
                raise ValueError(f'TPM row {line}: empty or duplicate gene ID')
            genes.add(gene)
            for name, index in indexes.items():
                value = float(row[index])
                if not math.isfinite(value) or value < 0:
                    raise ValueError(f'TPM row {line}, {name}: nonfinite or negative value')
                totals[name] += value
        if not genes:
            raise ValueError('TPM matrix has no gene rows')
        if not all(math.isfinite(value) for value in totals.values()):
            raise ValueError('TPM column sums overflowed')
    return len(genes), totals


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samplesheet', required=True, help='Local CSV path or public HTTP(S) CSV URL')
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--trimmer', choices=['fastp', 'trimgalore'], required=True)
    args = parser.parse_args()
    if not args.outdir.is_dir():
        parser.error(f'Output directory does not exist: {args.outdir}')
    rows = []

    def add(sample, check, status, detail):
        rows.append([sample, check, status, detail])

    def nonempty(sample, label, path):
        ok = path.is_file() and path.stat().st_size > 0
        add(sample, label, 'PASS' if ok else 'FAIL', str(path))
        return ok

    try:
        samples = load_samples(args.samplesheet)
        paired = sum(entry['paired'] for entry in samples.values())
        runs = sum(entry['runs'] for entry in samples.values())
        add('ALL', 'samplesheet', 'PASS',
            f'{runs} runs grouped into {len(samples)} samples; {paired} paired-end, {len(samples) - paired} single-end')
    except (OSError, ValueError, csv.Error) as error:
        samples = {}
        add('ALL', 'samplesheet', 'FAIL', str(error))

    if samples:
        tpm = args.outdir / 'tximport/salmon_merged.gene_tpm.tsv'
        try:
            genes, totals = check_tpm(tpm, samples)
            add('ALL', 'TPM_structure_and_values', 'PASS', f'{genes} genes; sample columns match')
            for name, total in totals.items():
                status = 'PASS' if abs(total - 1_000_000) <= 1 else 'WARN'
                add(name, 'TPM_sum', status, f'{total:.6f}; informational normalization check')
        except (OSError, ValueError, csv.Error) as error:
            add('ALL', 'TPM_structure_and_values', 'FAIL', str(error))

        for name, entry in samples.items():
            paired = entry['paired']
            counts = {}
            for stage, label in [('cat', 'merged_FASTQ'), (args.trimmer, 'trimmed_FASTQ')]:
                # Preserve compatibility with historical results that predate CAT_FASTQ.
                if stage == 'cat' and not (args.outdir / 'cat').is_dir():
                    add(name, label, 'WARN', 'cat/ absent; merged-read checks skipped for historical outputs')
                    continue
                paths = output_fastqs(args.outdir, name, paired, stage)
                try:
                    count = check_pair(*paths) if paired else check_single(paths[0])
                    counts[stage] = count
                    unit = 'pairs' if paired else 'reads'
                    detail = f'{count} {unit}; valid records'
                    if paired:
                        detail += ' and matching R1/R2 IDs'
                    if stage == 'cat':
                        detail += f"; {entry['runs']} input runs listed (input read totals not checked)"
                    add(name, label, 'PASS', detail)
                except (OSError, EOFError, ValueError, UnicodeError) as error:
                    add(name, label, 'FAIL', str(error))
            if 'cat' in counts and args.trimmer in counts:
                merged, trimmed = counts['cat'], counts[args.trimmer]
                add(name, 'trimmed_count_vs_merged', 'PASS' if trimmed <= merged else 'FAIL',
                    f'{trimmed} trimmed / {merged} merged; {trimmed / merged:.2%} retained')
            nonempty(name, 'Salmon_quant_sf', args.outdir / 'salmon' / name / 'quant.sf')
            directory = args.outdir / 'markduplicates'
            nonempty(name, 'marked_BAM', directory / f'{name}_marked.bam')
            indexes = [directory / f'{name}_marked.bai', directory / f'{name}_marked.bam.bai']
            index = next((p for p in indexes if p.is_file() and p.stat().st_size > 0), indexes[0])
            nonempty(name, 'marked_BAM_index', index)
            nonempty(name, 'duplication_metrics', directory / f'{name}_marked.metrics.txt')
        nonempty('ALL', 'MultiQC_report', args.outdir / 'multiqc/multiqc_report.html')

    report = args.outdir / 'validation_report.tsv'
    with report.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['sample', 'check', 'status', 'detail'])
        writer.writerows(rows)
    counts = {status: sum(row[2] == status for row in rows) for status in ('PASS', 'WARN', 'FAIL')}
    print(f"PASS={counts['PASS']} WARN={counts['WARN']} FAIL={counts['FAIL']}")
    print(f'Report: {report}')
    print('BAM, index, Salmon and report checks verify existence/nonempty files, not internal format integrity.')
    return 1 if counts['FAIL'] else 0


if __name__ == '__main__':
    raise SystemExit(main())
