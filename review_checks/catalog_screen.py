"""Screen a FULL downloaded GWAS Catalog association release against 37 SNPs and EUR proxies.

Input may be .tsv, .csv, .gz, or a .zip holding a single association table. A full
screen must contain at least 1,000,000 records; --allow-subset is for unit testing
against the saved 45 matched records and MUST NOT be presented as a full screen.
Reads large releases in chunks without querying the rate-limited API.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import zipfile
import pandas as pd

RS_PATTERN = re.compile(r'(?<![A-Za-z0-9])rs\d+(?!\d)', flags=re.IGNORECASE)

def hash_file(path: Path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(2**20), b''):
            h.update(b)
    return h.hexdigest()


def catalog_handle(path: Path, zip_member: str | None):
    if not path.is_file(): raise FileNotFoundError(path)
    if path.suffix.lower() == '.zip':
        z = zipfile.ZipFile(path)
        entries = [n for n in z.namelist() if n.lower().endswith(('.tsv', '.csv', '.txt', '.tsv.gz', '.csv.gz'))]
        if zip_member:
            if zip_member not in entries: raise ValueError(f'Archive lacks requested member {zip_member}')
            selected = zip_member
        elif len(entries) == 1:
            selected = entries[0]
        else:
            z.close(); raise ValueError(f'Archive has {len(entries)} tabular files: pass --zip-member')
        stream = z.open(selected)
        if selected.endswith('.gz'): stream = gzip.GzipFile(fileobj=stream)
        return stream,z,selected
    if path.name.endswith('.gz'): return gzip.open(path,'rb'),None,path.name
    return path.open('rb'),None,path.name


def run(catalog: Path, proxies: Path, output_dir: Path, allow_subset=False, zip_member=None, chunksize=50_000):
    proxy = pd.read_csv(proxies, dtype={'lookup_rsid':str,'instrument':str})
    needed = set(proxy.lookup_rsid.astype(str))
    lookup = {}
    for row in proxy.itertuples():
        lookup.setdefault(str(row.lookup_rsid),set()).add(str(row.instrument))
    output_dir.mkdir(parents=True,exist_ok=True)
    stream, z, member = catalog_handle(catalog, zip_member)
    try:
        line = stream.readline().decode('utf-8-sig')
        sep = '\t' if line.count('\t')>=line.count(',') else ','
    finally:
        stream.close()
        if z is not None: z.close()
    stream, z, member = catalog_handle(catalog, zip_member)
    rows, total = [],0
    try:
        read = pd.read_csv(stream,sep=sep,chunksize=chunksize,dtype=str,keep_default_na=False,low_memory=False)
        for chunk in read:
            if 'SNPS' not in chunk.columns:
                raise ValueError('Full GWAS Catalog file must contain SNPS column')
            total += len(chunk)
            for line in chunk.to_dict(orient='records'):
                matched = needed.intersection(set(RS_PATTERN.findall(line.get('SNPS','').lower())))
                if not matched: continue
                line['MATCHING_RSID'] = ';'.join(sorted(matched))
                line['MATCHING_INSTRUMENTS'] = ';'.join(sorted(set().union(*(lookup[rs] for rs in matched))))
                rows.append(line)
    finally:
        stream.close()
        if z is not None:z.close()
    if total < 1_000_000 and not allow_subset:
        raise RuntimeError(f'Input contains only {total:,} records; NOT a full Catalog release. Use --allow-subset for tests.')
    out = pd.DataFrame(rows)
    out.to_csv(output_dir/'catalog_screen_matches.tsv',sep='\t',index=False)
    report = {'catalog_path':str(catalog),'catalog_member':member,'catalog_sha256':hash_file(catalog), 'association_records_screened':total,'full_release_screen':bool(total>=1_000_000), 'lookup_rsid_count':len(needed),'matched_association_records':len(out),'matched_unique_rsids':len(set().union(*(set(x.split(';')) for x in out.MATCHING_RSID)) if len(out) else set()),'output':'catalog_screen_matches.tsv'}
    (output_dir/'catalog_screen_report.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    if not report['full_release_screen']:
        print('UNIT-TEST SUBSET: Do not report this as an exhaustive screen.')

if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    base=Path(__file__).resolve().parent
    ap.add_argument('--catalog',type=Path,required=True)
    ap.add_argument('--proxies',type=Path,default=base/'eur_proxies_r2_0.8.csv')
    ap.add_argument('--out-dir',type=Path,default=base)
    ap.add_argument('--allow-subset',action='store_true')
    ap.add_argument('--zip-member',default=None)
    args=ap.parse_args()
    run(args.catalog.resolve(),args.proxies.resolve(),args.out_dir.resolve(),args.allow_subset,args.zip_member)
