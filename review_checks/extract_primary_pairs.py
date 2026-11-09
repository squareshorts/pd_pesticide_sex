"""Recreate the four 37-SNP input slices for reciprocal_sex.py from the ORIGINAL full GWAS files.

Usage from pd_pesticide_sex/review_checks:
  python extract_primary_pairs.py --data-dir ../data

Never selects variants based on outcome associations. Matches the previously selected 37
instruments by chromosome, position, and allele pair, checks orientation, and produces
an audit of missing / nonunique matches. Source files are not modified.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
import pandas as pd

COMPLEMENT = str.maketrans('ACGT', 'TGCA')
SEXES = ('female', 'male')
LEADS = ('rs62092069', 'rs12456492', 'rs12955331', 'rs1944418', 'rs8087199')
EXPOSURE_FIELDS = ['variant', 'beta', 'se', 'pval', 'AC', 'n_complete_samples', 'low_confidence_variant']
OUTCOME_FIELDS = ['ID', 'Allele1', 'Allele2', 'Freq1', 'Effect', 'StdErr', 'P-value']


def parse_variant(variant: str):
    parts = str(variant).split(':')
    if len(parts) != 4:
        raise ValueError(f'Unexpected Neale variant syntax: {variant!r}')
    chrom, pos, ref, alt = parts
    return chrom.removeprefix('chr'), int(pos), ref.upper(), alt.upper()


def align_variant(source_variant: str, target_variant: str):
    """Return +1/-1 for source ALT effect relative to target ALT; None for mismatch.

    Exact allele order is preferred. Reverse-complement pairs are accepted only when
    strand orientation is unambiguous. A/T and C/G swapped pairs are ambiguous and
    are not silently aligned.
    """
    c, p, ref, alt = parse_variant(source_variant)
    d, q, target_ref, target_alt = parse_variant(target_variant)
    if (c, p) != (d, q):
        return None
    if (ref, alt) == (target_ref, target_alt):
        return 1
    palindromic = {ref, alt} in ({'A', 'T'}, {'C', 'G'})
    if (ref, alt) == (target_alt, target_ref):
        return None if palindromic else -1
    cr, ca = ref.translate(COMPLEMENT), alt.translate(COMPLEMENT)
    if (cr, ca) == (target_ref, target_alt):
        return None if palindromic else 1
    if (cr, ca) == (target_alt, target_ref):
        return None if palindromic else -1
    return None


def write_table(frame: pd.DataFrame, path: Path, columns: list[str]):
    frame.reindex(columns=columns).to_csv(path, sep='\t', index=False)
    print(f'{path.name}: {len(frame)} variants')


def run(panel: Path, data_dir: Path, out: Path, chunksize: int = 250_000):
    if not panel.is_file():
        raise FileNotFoundError(panel)
    out.mkdir(parents=True, exist_ok=True)
    p = pd.read_csv(panel, dtype={'chrom': str})
    assert len(p) == 37 and p.rsid.is_unique and p.variant.is_unique, 'Expected frozen 37-SNP panel'
    target = {tuple(parse_variant(v)[:2]): str(v) for v in p.variant}
    assert len(target) == 37, 'Nonunique genomic positions in selected panel'
    audit = {'panel': str(panel), 'source_files': {}, 'exposure': {}, 'outcome': {}, 'leads': list(LEADS)}

    for sex in SEXES:
        full = data_dir / f'22614_1.gwas.imputed_v3.{sex}.tsv.bgz'
        if not full.is_file():
            raise FileNotFoundError(f'Missing full Neale GWAS: {full}')
        audit['source_files'][str(full.name)] = int(full.stat().st_size)
        hits = {}
        total = 0
        for chunk in pd.read_csv(full, sep='\t', compression='gzip', chunksize=chunksize, usecols=EXPOSURE_FIELDS):
            total += len(chunk)
            for rec in chunk.itertuples(index=False):
                record = rec._asdict()
                try:
                    chrom, pos, _, _ = parse_variant(record['variant'])
                except ValueError:
                    continue
                key = (chrom, pos)
                if key not in target:
                    continue
                hits.setdefault(key, []).append(record)
        rows = []
        missing = []
        for key, canonical in target.items():
            matches = []
            for record in hits.get(key, []):
                orientation = align_variant(record['variant'], canonical)
                if orientation is not None:
                    matches.append((record, orientation))
            if len(matches) != 1:
                missing.append({'variant': canonical, 'reason': 'absent_or_ambiguous', 'position_hits': len(hits.get(key, [])), 'allele_matches': len(matches)})
                continue
            record, orientation = matches[0]
            n = float(record['n_complete_samples'])
            ac = float(record['AC'])
            record['beta'] = float(record['beta']) * orientation
            record['AC'] = ac if orientation == 1 else 2 * n - ac
            record['variant'] = canonical
            rows.append(record)
        if missing:
            audit['exposure'][sex] = {'rows_scanned': total, 'missing_or_ambiguous': missing}
            (out / 'extraction_audit.json').write_text(json.dumps(audit, indent=2))
            raise RuntimeError(f'Unable to recover all 37 exposure records for {sex}; see extraction_audit.json')
        result = pd.DataFrame(rows)
        write_table(result, out / f'neale_22614_1_{sex}_37.tsv', EXPOSURE_FIELDS)
        audit['exposure'][sex] = {'rows_scanned': total, 'matched_37': len(result), 'n_reoriented': sum(align_variant(hits[k][0]['variant'], target[k]) == -1 for k in target)}

    wanted = set(p.rsid) | set(LEADS)
    for sex in SEXES:
        full = data_dir / f'{sex.upper()}_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt.gz'
        if not full.is_file():
            raise FileNotFoundError(f'Missing full PD outcome GWAS: {full}')
        audit['source_files'][str(full.name)] = int(full.stat().st_size)
        found = []
        scanned = 0
        for chunk in pd.read_csv(full, sep='\t', compression='gzip', chunksize=chunksize, usecols=OUTCOME_FIELDS):
            scanned += len(chunk)
            chunk = chunk[chunk.ID.astype(str).isin(wanted)]
            if not chunk.empty:
                found.append(chunk)
        fullhits = pd.concat(found, ignore_index=True) if found else pd.DataFrame(columns=OUTCOME_FIELDS)
        if fullhits.ID.duplicated().any():
            raise RuntimeError(f'{sex} PD GWAS has duplicate requested IDs')
        found_main = fullhits[fullhits.ID.isin(set(p.rsid))]
        found_leads = fullhits[fullhits.ID.isin(set(LEADS))]
        write_table(found_main, out / f'ipdgc_{sex.upper()}_37.tsv', OUTCOME_FIELDS)
        write_table(found_leads, out / f'ipdgc_{sex.upper()}_RIT2.tsv', OUTCOME_FIELDS)
        audit['outcome'][sex] = {'rows_scanned': scanned, 'n_primary_found': len(found_main), 'missing_primary': sorted(set(p.rsid) - set(found_main.ID)), 'n_rit2_found': len(found_leads), 'missing_rit2': sorted(set(LEADS) - set(found_leads.ID))}
    (out / 'extraction_audit.json').write_text(json.dumps(audit, indent=2))
    print('Extraction complete. Run reciprocal_sex.py from this directory.')


if __name__ == '__main__':
    base = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data-dir', type=Path, default=base.parent / 'data')
    ap.add_argument('--panel', type=Path, default=base / 'primary_panel_37.csv')
    ap.add_argument('--out-dir', type=Path, default=base)
    ap.add_argument('--chunksize', type=int, default=250_000)
    args = ap.parse_args()
    run(args.panel.resolve(), args.data_dir.resolve(), args.out_dir.resolve(), args.chunksize)
