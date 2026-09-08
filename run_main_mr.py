import csv,gzip,math,os,json
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

BASE='/mnt/data'
settings={
 'female': {
  'inst': f'{BASE}/main_instr/female_instruments.csv',
  'out': f'{BASE}/FEMALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt.gz',
 },
 'male': {
  'inst': f'{BASE}/main_instr/male_instruments.csv',
  'out': f'{BASE}/MALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt.gz',
 }
}

def complement(a):
    tr=str.maketrans('ATCG','TAGC')
    return a.translate(tr) if all(c in 'ATCG' for c in a) else None

def is_pal(a,b):
    return len(a)==1 and len(b)==1 and {a,b} in ({'A','T'},{'C','G'})

def extract_outcome(path, ids):
    wanted=set(ids); found={}
    with gzip.open(path,'rt',newline='') as f:
        rd=csv.DictReader(f, delimiter='\t')
        for r in rd:
            rs=r.get('ID')
            if rs in wanted:
                found[rs]=r
                if len(found)==len(wanted): break
    return found

def harmonize(inst, outmap):
    rows=[]; drops=[]
    for _,x in inst.iterrows():
        rs=x.rsid
        if rs not in outmap:
            drops.append((rs,'absent_outcome')); continue
        r=outmap[rs]
        ea=str(x.effect_allele).upper(); oa=str(x.other_allele).upper()
        yea=str(r['Allele1']).upper(); yoa=str(r['Allele2']).upper()
        try:
            by=float(r['Effect']); sy=float(r['StdErr']); py=float(r['P-value']); fy=float(r['Freq1'])
        except:
            drops.append((rs,'bad_outcome_numeric')); continue
        bx=float(x.beta); sx=float(x.se); ex=float(x.eaf)
        if not all(math.isfinite(v) for v in [by,sy,fy,bx,sx,ex]) or sy<=0 or sx<=0:
            drops.append((rs,'nonfinite')); continue
        flip=None; strand='same'
        # exact allele orientation
        if ea==yea and oa==yoa: flip=False
        elif ea==yoa and oa==yea: flip=True
        else:
            cea=complement(ea); coa=complement(oa)
            if cea is not None and coa is not None:
                if cea==yea and coa==yoa: flip=False; strand='complement'
                elif cea==yoa and coa==yea: flip=True; strand='complement'
        if flip is None:
            drops.append((rs,f'alleles_mismatch:{ea}/{oa}:{yea}/{yoa}')); continue
        # Palindromes: discard if exposure MAF ambiguous. Otherwise frequency-check the chosen orientation.
        if is_pal(ea,oa):
            if 0.42 < ex < 0.58:
                drops.append((rs,'palindromic_ambiguous_maf')); continue
            aligned_fy = 1-fy if flip else fy
            # If chosen orientation gives incompatible allele freq but reverse gives compatible, reverse.
            d1=abs(ex-aligned_fy)
            alt_fy=fy if flip else 1-fy
            d2=abs(ex-alt_fy)
            if d1>0.20 and d2<d1:
                flip=not flip; aligned_fy=alt_fy; d1=d2
            if d1>0.20:
                drops.append((rs,f'palindromic_freq_mismatch:{d1:.3f}')); continue
        by_al=-by if flip else by
        fy_al=1-fy if flip else fy
        rows.append({**x.to_dict(),'beta_outcome':by_al,'se_outcome':sy,'p_outcome':py,
                     'eaf_outcome_aligned':fy_al,'outcome_effect_allele_raw':yea,'outcome_other_allele_raw':yoa,
                     'flip_outcome':bool(flip),'strand_action':strand})
    return pd.DataFrame(rows), drops

def ivw(df, random=True):
    bx=df.beta.to_numpy(float); by=df.beta_outcome.to_numpy(float); sy=df.se_outcome.to_numpy(float)
    w=1/sy**2
    den=np.sum(w*bx*bx); b=np.sum(w*bx*by)/den
    se_fixed=np.sqrt(1/den)
    resid=by-b*bx; Q=np.sum(w*resid**2); dof=len(df)-1; q_p=stats.chi2.sf(Q,dof) if dof>0 else np.nan
    phi=max(1.0,Q/dof) if random and dof>0 else 1.0
    se=se_fixed*np.sqrt(phi)
    z=b/se; p=2*stats.norm.sf(abs(z))
    return dict(beta=b,se=se,p=p,ci_low=b-1.96*se,ci_high=b+1.96*se,OR=np.exp(b),OR_low=np.exp(b-1.96*se),OR_high=np.exp(b+1.96*se),Q=Q,Q_df=dof,Q_p=q_p,phi=phi,se_fixed=se_fixed)

def egger(df):
    bx=df.beta.to_numpy(float); by=df.beta_outcome.to_numpy(float); sy=df.se_outcome.to_numpy(float)
    X=sm.add_constant(bx); w=1/sy**2
    fit=sm.WLS(by,X,weights=w).fit()
    # Multiplicative random-effects style: WLS covariance already scaled by residual SSR/df.
    return dict(intercept=float(fit.params[0]), intercept_se=float(fit.bse[0]), intercept_p=float(fit.pvalues[0]),
                beta=float(fit.params[1]),se=float(fit.bse[1]),p=float(fit.pvalues[1]),
                ci_low=float(fit.conf_int()[1,0]),ci_high=float(fit.conf_int()[1,1]),
                OR=float(np.exp(fit.params[1])),OR_low=float(np.exp(fit.conf_int()[1,0])),OR_high=float(np.exp(fit.conf_int()[1,1])))

def weighted_median_point(df, by_override=None):
    bx=df.beta.to_numpy(float); by=df.beta_outcome.to_numpy(float) if by_override is None else by_override; sy=df.se_outcome.to_numpy(float)
    ratio=by/bx; weights=(bx**2)/(sy**2)
    idx=np.argsort(ratio); r=ratio[idx]; w=weights[idx]; cw=np.cumsum(w)/np.sum(w)
    j=np.searchsorted(cw,0.5)
    if j==0: return float(r[0])
    # linear interpolate around median cumulative probability
    c0,c1=cw[j-1],cw[j]; r0,r1=r[j-1],r[j]
    if c1==c0: return float(r[j])
    return float(r0+(0.5-c0)/(c1-c0)*(r1-r0))

def weighted_median(df, B=10000, seed=20260908):
    b=weighted_median_point(df)
    rng=np.random.default_rng(seed)
    by=df.beta_outcome.to_numpy(float); sy=df.se_outcome.to_numpy(float)
    sims=np.empty(B)
    for i in range(B): sims[i]=weighted_median_point(df, rng.normal(by,sy))
    se=float(np.std(sims,ddof=1)); p=2*stats.norm.sf(abs(b/se)) if se>0 else np.nan
    lo,hi=np.quantile(sims,[.025,.975])
    return dict(beta=b,se=se,p=p,ci_low=float(lo),ci_high=float(hi),OR=np.exp(b),OR_low=np.exp(lo),OR_high=np.exp(hi))

def robust_ivw_huber(df):
    bx=df.beta.to_numpy(float); by=df.beta_outcome.to_numpy(float); sy=df.se_outcome.to_numpy(float)
    # standardize by outcome SE; robust regression through origin
    X=(bx/sy)[:,None]; Y=by/sy
    fit=sm.RLM(Y,X,M=sm.robust.norms.HuberT()).fit()
    b=float(fit.params[0]); se=float(fit.bse[0]); p=2*stats.norm.sf(abs(b/se));
    return dict(beta=b,se=se,p=p,ci_low=b-1.96*se,ci_high=b+1.96*se,OR=np.exp(b),OR_low=np.exp(b-1.96*se),OR_high=np.exp(b+1.96*se))

def loo(df):
    rows=[]
    for i,r in df.iterrows():
        sub=df.drop(index=i)
        z=ivw(sub,True)
        rows.append({'left_out':r.rsid,**{k:z[k] for k in ['beta','se','p','ci_low','ci_high']}})
    return pd.DataFrame(rows)

allres={}; harmonized={}
for sex,s in settings.items():
    inst=pd.read_csv(s['inst'])
    print(sex,'extracting',len(inst),'outcome SNPs...')
    out=extract_outcome(s['out'], inst.rsid)
    print(sex,'found',len(out),'of',len(inst))
    dat,drops=harmonize(inst,out)
    print(sex,'harmonized',len(dat),'drops',drops)
    dat.to_csv(f'{BASE}/{sex}_harmonized_main.csv',index=False)
    lr=loo(dat); lr.to_csv(f'{BASE}/{sex}_loo_main.csv',index=False)
    allres[sex]={
      'n_initial':len(inst),'n_outcome_found':len(out),'n_harmonized':len(dat),'drops':drops,
      'mean_F':float(dat.F.mean()),'min_F':float(dat.F.min()),'max_F':float(dat.F.max()),
      'ivw_mre':ivw(dat,True),'ivw_fixed':ivw(dat,False),'egger':egger(dat),'weighted_median':weighted_median(dat),'robust_huber':robust_ivw_huber(dat),
      'loo_beta_min':float(lr.beta.min()),'loo_beta_max':float(lr.beta.max()),'loo_p_min':float(lr.p.min()),'loo_p_max':float(lr.p.max())
    }
    harmonized[sex]=dat

f=allres['female']['ivw_mre']; m=allres['male']['ivw_mre']
diff=f['beta']-m['beta']; sed=math.sqrt(f['se']**2+m['se']**2); z=diff/sed; p=2*stats.norm.sf(abs(z))
interaction={'beta_difference_female_minus_male':diff,'se_difference':sed,'z':z,'p':p,'ci_low':diff-1.96*sed,'ci_high':diff+1.96*sed}
# Fixed effect meta of sex-specific IVWs for descriptive common estimate
w_f=1/f['se']**2; w_m=1/m['se']**2; bmeta=(w_f*f['beta']+w_m*m['beta'])/(w_f+w_m); semeta=math.sqrt(1/(w_f+w_m)); pmeta=2*stats.norm.sf(abs(bmeta/semeta))
meta={'beta':bmeta,'se':semeta,'p':pmeta,'ci_low':bmeta-1.96*semeta,'ci_high':bmeta+1.96*semeta,'OR':math.exp(bmeta),'OR_low':math.exp(bmeta-1.96*semeta),'OR_high':math.exp(bmeta+1.96*semeta)}
summary={'sex_specific':allres,'sex_interaction':interaction,'fixed_effect_meta_of_sex_estimates':meta}
with open(f'{BASE}/main_mr_summary.json','w') as o: json.dump(summary,o,indent=2)
print(json.dumps(summary,indent=2))
