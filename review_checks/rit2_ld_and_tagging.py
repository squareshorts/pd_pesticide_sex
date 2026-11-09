"""Recreate RIT2 LD and regional PD-tagging calculations using 1000G EUR genotypes.

Requires data/EUR.bed,bim,fam plus the *_RIT2.tsv outcome extracts produced by
extract_primary_pairs.py. Run from the review_checks directory. Genotype dosages
refer to the A1 allele in EUR.bim. Output coefficients are explicitly reoriented.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import ld_tools as L

TARGET = 'rs62092069'
LEADS = ['rs12456492', 'rs12955331', 'rs1944418', 'rs8087199']
COMP = str.maketrans('ACGT', 'TGCA')


def align_effect(effect: float, a: str, b: str, to_a: str, to_b: str):
    a, b, to_a, to_b = (str(x).upper() for x in (a, b, to_a, to_b))
    if (a, b) == (to_a, to_b): return float(effect)
    if (a, b) == (to_b, to_a): return -float(effect)
    if (a.translate(COMP), b.translate(COMP)) == (to_a, to_b): return float(effect)
    if (a.translate(COMP), b.translate(COMP)) == (to_b, to_a): return -float(effect)
    raise ValueError(f'Unable to harmonize [{a}/{b}] with [{to_a}/{to_b}]')


def run():
    base = Path(__file__).resolve().parent
    panel = pd.read_csv(base / 'primary_panel_37.csv')
    target = panel.loc[panel.rsid == TARGET]
    if len(target) != 1: raise ValueError('Missing target instrument in frozen panel')
    effect_allele = str(target.iloc[0].effect_allele)
    wanted = [TARGET] + LEADS
    bim = L.bim_lookup(wanted)
    by_id = {}
    for item in bim:
        by_id.setdefault(item[2], []).append(item)
    for rs in wanted:
        if len(by_id.get(rs, [])) != 1:
            raise ValueError(f'Expected exactly one reference genotype for {rs}, got {len(by_id.get(rs, []))}')
    variants = {rs: by_id[rs][0] for rs in wanted}
    G = L.genotypes([variants[rs][0] for rs in wanted])
    gt = G[variants[TARGET][0]]
    ok = np.isfinite(gt)
    if ok.sum() < 20: raise RuntimeError('Insufficient nonmissing genotypes at target')
    result_ld = {}
    result_tag = {}
    for rs in LEADS:
        item = variants[rs]
        lead = G[item[0]]
        mask = np.isfinite(gt) & np.isfinite(lead)
        if mask.sum() < 20: raise RuntimeError(f'Insufficient paired genotypes for {rs}')
        sd_target = float(np.std(gt[mask], ddof=0))
        sd_lead = float(np.std(lead[mask], ddof=0))
        r_a1 = float(np.corrcoef(gt[mask], lead[mask])[0, 1])
        if not np.isfinite(r_a1) or sd_target == 0: raise RuntimeError(f'Undefined LD for {rs}')
        result_ld[rs] = {'pos': int(item[3]), 'r2': r_a1**2, 'signed_r_A1_A1': r_a1, 'n_pairs': int(mask.sum())}
        target_ref = variants[TARGET]
        # Convert correlation of target reference A1 dosage to panel effect-allele dosage.
        panel_sign = align_effect(1.0, target_ref[4], target_ref[5], effect_allele, target_ref[5] if target_ref[4] == effect_allele else target_ref[4])
        for sex in ('FEMALE', 'MALE'):
            file = base / f'ipdgc_{sex}_RIT2.tsv'
            if not file.is_file(): raise FileNotFoundError(f'Missing {file}; run extract_primary_pairs.py')
            outcome = pd.read_csv(file, sep='\t').set_index('ID')
            if TARGET not in outcome.index or rs not in outcome.index:
                result_tag[f'{sex}_{rs}'] = {'status': 'missing_outcome'}
                continue
            t, lead_out = outcome.loc[TARGET], outcome.loc[rs]
            observed_target = align_effect(t['Effect'], t['Allele1'], t['Allele2'], effect_allele, target_ref[5] if effect_allele == target_ref[4] else target_ref[4])
            beta_lead_refa1 = align_effect(lead_out['Effect'], lead_out['Allele1'], lead_out['Allele2'], item[4], item[5])
            predicted = panel_sign * r_a1 * beta_lead_refa1 * sd_lead / sd_target
            result_tag[f'{sex}_{rs}'] = {'observed_beta_exposure_allele': observed_target, 'observed_se': float(t['StdErr']), 'observed_p': float(t['P-value']), 'lead_beta_reference_A1': beta_lead_refa1, 'lead_p': float(lead_out['P-value']), 'signed_r_target_effect_vs_lead_A1': panel_sign*r_a1, 'expected_beta_from_LD': predicted, 'r2': r_a1**2}
    (base/'rit2_ld_reproduced.json').write_text(json.dumps(result_ld, indent=2))
    (base/'rit2_tagging_reproduced.json').write_text(json.dumps(result_tag, indent=2))
    print('RIT2 LD and tagging calculations complete; see rit2_*_reproduced.json')

if __name__ == '__main__':
    run()
