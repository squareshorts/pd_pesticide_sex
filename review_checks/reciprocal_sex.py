"""Reviewer 4.5 / Reviewer 2: same-SNP sex comparison, cross-sex replication, and reciprocal-selection MR
for UKB 22614_1 ("worked with pesticides: sometimes") -> PD (IPDGC, no UKB).
Inputs: the 37 primary instruments (results/22614_1_{sex}_harmonized.csv), their rows extracted from BOTH full
sex-specific Neale files (extract step) and from BOTH IPDGC outcome files; LD from data/EUR.*.
Outcome harmonisation reproduces scripts/prepare.py harmonize() exactly."""
import json, numpy as np, pandas as pd
from scipy.stats import norm, binomtest, linregress
import ld_tools as L
comp=str.maketrans('ACGT','TGCA'); SEX=['female','male']
p=pd.read_csv('primary_panel_37.csv')
nx={s:pd.read_csv(f'neale_22614_1_{s}_37.tsv',sep='\t').set_index('variant') for s in SEX}
oy={s:pd.read_csv(f'ipdgc_{s.upper()}_37.tsv',sep='\t').set_index('ID') for s in SEX}

def harmonize_outcome(a,b,eaf,y):
    c,e=str(y['Allele1']).upper(),str(y['Allele2']).upper(); fy,by,sy=float(y['Freq1']),float(y['Effect']),float(y['StdErr'])
    pal={a,b} in [{'A','T'},{'C','G'}]; flip=None
    for ac,bc,st in [(a,b,False),(a.translate(comp),b.translate(comp),True)]:
        if (ac,bc)==(c,e): flip=False; break
        if (ac,bc)==(e,c): flip=True; break
    if flip is None: return None,'allele_mismatch'
    if not np.isfinite([fy,by,sy]).all() or sy<=0 or not 0<fy<1: return None,'invalid_outcome'
    if pal:
        if .42<=eaf<=.58 or .42<=fy<=.58: return None,'palindromic_ambiguous_frequency'
        af=1-fy if flip else fy
        if (eaf<.5)!=(af<.5): flip=not flip
        if abs(eaf-(1-fy if flip else fy))>.20: return None,'palindromic_frequency_mismatch'
    return dict(beta_PD=-by if flip else by,se_PD=sy,p_PD=float(y['P-value'])),None

rows=[];excl=[]
for _,r in p.iterrows():
    for s in SEX:
        base=dict(rsid=r.rsid,origin=r.origin,sex=s,variant=r.variant,effect_allele=r.effect_allele,other_allele=r.other_allele)
        if r.variant not in nx[s].index: excl.append({**base,'reason':'absent_exposure'}); continue
        x=nx[s].loc[r.variant]   # Neale beta refers to the ALT allele = effect_allele in the primary panel
        if not (np.isfinite(x.beta) and np.isfinite(x.se) and x.se>0): excl.append({**base,'reason':'invalid_exposure'}); continue
        eaf=x.AC/(2*x.n_complete_samples)
        if r.rsid not in oy[s].index: excl.append({**base,'reason':'absent_outcome'}); continue
        h,why=harmonize_outcome(r.effect_allele,r.other_allele,eaf,oy[s].loc[r.rsid])
        if why: excl.append({**base,'reason':why}); continue
        rows.append({**base,'beta_X':x.beta,'se_X':x.se,'p_X':x.pval,'eaf':eaf,'low_confidence':str(x.low_confidence_variant),'n':x.n_complete_samples,**h})
d=pd.DataFrame(rows); pd.DataFrame(excl).to_csv('reciprocal_exclusions.csv',index=False)
# sanity: own-sex rows must reproduce the archived harmonised panel
own=d[d.sex==d.origin].merge(p,on='rsid'); chk=dict(max_abs_dbeta_X=float((own.beta_X-own.beta).abs().max()),max_abs_dbeta_PD=float((own.beta_PD-own.beta_outcome).abs().max()))
paired=d.groupby('rsid').sex.nunique(); both=paired[paired==2].index; d=d[d.rsid.isin(both)].copy()

# union LD clump (r2>0.001 within 10 Mb), priority = smallest discovery p (own-sex selection p)
pri=d[d.sex==d.origin][['rsid','p_X']].sort_values('p_X')
bim=L.bim_lookup(pri.rsid); info={r[2]:r for r in bim}; G=L.genotypes([r[0] for r in bim])
kept=[];pruned={}
for rs in pri.rsid:
    i,c,_,pos,_,_=info[rs]; clash=[k for k in kept if info[k][1]==c and abs(info[k][3]-pos)<=10_000_000 and L.r2(G[i],G[info[k][0]])>0.001]
    if clash: pruned[rs]=clash
    else: kept.append(rs)
d.to_csv('reciprocal_all_pairs.csv',index=False)
cp=d[d.rsid.isin(kept)].copy(); cp.to_csv('reciprocal_common_panel.csv',index=False)

def ivw(g):
    x,y,s=g.beta_X.values,g.beta_PD.values,g.se_PD.values; w=1/s**2; den=(w*x*x).sum(); b=(w*x*y).sum()/den
    Q=(w*(y-b*x)**2).sum(); phi=max(1,Q/(len(g)-1)); se=np.sqrt(phi/den); z=b/se
    F=(x/g.se_X.values)**2
    return dict(K=len(g),beta=b,se=se,ci95=[b-1.959964*se,b+1.959964*se],p=2*norm.sf(abs(z)),Q=Q,Q_p=float(__import__('scipy.stats',fromlist=['chi2']).chi2.sf(Q,len(g)-1)),F_mean=float(F.mean()),F_median=float(np.median(F)),F_min=float(F.min()),n_F_gt_10=int((F>10).sum()))
def contrast(a,b):
    dl=a['beta']-b['beta']; se=np.hypot(a['se'],b['se']); return dict(delta=dl,se=se,ci95=[dl-1.959964*se,dl+1.959964*se],p=2*norm.sf(abs(dl/se)))
res={'sanity_own_sex_vs_archived':chk,'n_pairs':int(len(both)),'n_after_union_clump':len(kept),'pruned_by_union_clump':pruned}
# 1. same-SNP panel
f=ivw(cp[cp.sex=='female']); m=ivw(cp[cp.sex=='male']); res['common_panel']={'female':f,'male':m,'female_minus_male':contrast(f,m)}
# 2. reciprocal selection: each sex's MR uses only instruments discovered in the OTHER sex (exposure effects unaffected by own-sample winner's curse)
rf=ivw(cp[(cp.sex=='female')&(cp.origin=='male')]); rm=ivw(cp[(cp.sex=='male')&(cp.origin=='female')])
res['reciprocal_selection']={'female_using_male_selected':rf,'male_using_female_selected':rm,'female_minus_male':contrast(rf,rm)}
# 3. cross-sex replication of SNP-exposure associations
rep={}
for o in SEX:
    other=[s for s in SEX if s!=o][0]
    a=cp[(cp.origin==o)&(cp.sex==o)].set_index('rsid'); b=cp[(cp.origin==o)&(cp.sex==other)].set_index('rsid').loc[a.index]
    conc=int((np.sign(a.beta_X)==np.sign(b.beta_X)).sum()); n=len(a)
    lr=linregress(a.beta_X,b.beta_X)
    # inverse-variance weighted slope of replication effects on discovery effects
    w=1/b.se_X**2; slope_w=float((w*a.beta_X*b.beta_X).sum()/(w*a.beta_X**2).sum())
    zrep=b.beta_X/b.se_X*np.sign(a.beta_X)  # replication z oriented to discovery direction
    rep[f'{o}_selected_tested_in_{other}']=dict(N=n,sign_concordant=conc,binomial_p=binomtest(conc,n,.5).pvalue,
        n_same_direction_p_lt_05=int(((zrep>1.959964)).sum()),expected_under_null=round(0.025*n,2),
        mean_oriented_z=float(zrep.mean()),sum_z_over_sqrtN=float(zrep.sum()/np.sqrt(n)),p_sum_z_one_sided=float(norm.sf(zrep.sum()/np.sqrt(n))),
        ols_slope=lr.slope,ols_r=lr.rvalue,ols_p=lr.pvalue,ivw_slope=slope_w,
        other_sex_F_mean=float(((b.beta_X/b.se_X)**2).mean()))
res['replication']=rep
def clean(o):
    if isinstance(o,dict): return {k:clean(v) for k,v in o.items()}
    if isinstance(o,(list,tuple)): return [clean(v) for v in o]
    if isinstance(o,(np.floating,float)): return float(o)
    if isinstance(o,(np.integer,)): return int(o)
    return o
json.dump(clean(res),open('reciprocal_summary.json','w'),indent=2); print(json.dumps(clean(res),indent=1))
