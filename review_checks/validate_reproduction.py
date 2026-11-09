"""Compare outputs regenerated from original public inputs with submitted compact results.

Run after extract_primary_pairs.py, reciprocal_sex.py, rit2_ld_and_tagging.py,
proxies.py, and catalog_screen.py. Fails on a material discrepancy. The Catalog
check only passes when a full >=1M-record release was actually screened.
"""
from pathlib import Path
import json
import math
import pandas as pd

BASE = Path(__file__).resolve().parent

def close(a,b,tol):
    if not (math.isfinite(float(a)) and math.isfinite(float(b)) and abs(float(a)-float(b))<=tol):
        raise AssertionError(f'{a} != {b} within {tol}')


def main():
    summary=json.loads((BASE/'reciprocal_summary.json').read_text())
    if summary['n_pairs']!=35 or summary['n_after_union_clump']!=32:raise AssertionError('Reciprocal panel size changed')
    close(summary['common_panel']['female']['beta'],-1.20457191224,2e-6)
    close(summary['common_panel']['female']['se'],1.80891866799,2e-6)
    close(summary['common_panel']['male']['beta'],-0.06610335840,2e-6)
    close(summary['common_panel']['male']['se'],0.82273098226,2e-6)
    r=summary['replication']
    if r['female_selected_tested_in_male']['sign_concordant']+r['male_selected_tested_in_female']['sign_concordant']!=17:
        raise AssertionError('Cross-sex concordance differs from submitted results')
    ld=json.loads((BASE/'rit2_ld_reproduced.json').read_text())
    old=json.loads((BASE/'rit2_ld_results.json').read_text())
    for rs,rec in old.items():close(ld[rs]['r2'],rec['r2'],0.015)
    tag=json.loads((BASE/'rit2_tagging_reproduced.json').read_text())
    oldtag=json.loads((BASE/'rit2_tagging_expectation.json').read_text())
    for name,rec in oldtag.items():
        if tag.get(name,{}).get('status')=='missing_outcome':raise AssertionError(f'Missing outcome for {name}')
        close(tag[name]['expected_beta_from_LD'],rec['expected_beta_from_LD'],0.015)
    proxy=pd.read_csv(BASE/'eur_proxies_r2_0.8.csv')
    if proxy.lookup_rsid.nunique()!=274:raise AssertionError('Proxy count differs')
    catalog=json.loads((BASE/'catalog_screen_report.json').read_text())
    if not catalog['full_release_screen']:raise AssertionError('Catalog was screened only on a subset')
    if catalog['matched_association_records'] !=45:raise AssertionError('Catalog match count differs')
    print('PASS: common panel, sign concordance, RIT2 LD, proxy count, complete Catalog screen')

if __name__=='__main__':main()
