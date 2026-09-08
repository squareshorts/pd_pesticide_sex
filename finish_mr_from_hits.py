import pandas as pd, numpy as np, math, json
from scipy import stats
import statsmodels.api as sm
B='/mnt/data'

def comp(a):
 tr=str.maketrans('ATCG','TAGC'); return a.translate(tr) if all(c in 'ATCG' for c in a) else None

def pal(a,b): return len(a)==len(b)==1 and {a,b} in ({'A','T'},{'C','G'})

def harmonize(inst,out):
 om={r.ID:r for _,r in out.iterrows()}; rows=[]; drops=[]
 for _,x in inst.iterrows():
  if x.rsid not in om: drops.append([x.rsid,'absent_outcome']); continue
  r=om[x.rsid]; ea=str(x.effect_allele).upper(); oa=str(x.other_allele).upper(); yea=str(r.Allele1).upper(); yoa=str(r.Allele2).upper(); ex=float(x.eaf); fy=float(r.Freq1)
  flip=None; strand='same'
  if ea==yea and oa==yoa: flip=False
  elif ea==yoa and oa==yea: flip=True
  else:
   ce,co=comp(ea),comp(oa)
   if ce and co:
    if ce==yea and co==yoa: flip=False; strand='complement'
    elif ce==yoa and co==yea: flip=True; strand='complement'
  if flip is None: drops.append([x.rsid,f'mismatch:{ea}/{oa}:{yea}/{yoa}']); continue
  if pal(ea,oa):
   if .42<ex<.58: drops.append([x.rsid,'palindromic_ambiguous_maf']); continue
   aligned=1-fy if flip else fy; d1=abs(ex-aligned); alt=fy if flip else 1-fy; d2=abs(ex-alt)
   if d1>.2 and d2<d1: flip=not flip; aligned=alt; d1=d2
   if d1>.2: drops.append([x.rsid,f'pal_freq_mismatch:{d1:.3f}']); continue
  by=float(r.Effect); sy=float(r.StdErr)
  rows.append({**x.to_dict(),'beta_outcome':-by if flip else by,'se_outcome':sy,'eaf_outcome_aligned':1-fy if flip else fy,'flip_outcome':flip,'strand_action':strand})
 return pd.DataFrame(rows),drops

def ivw(d,random=True):
 bx=d.beta.values.astype(float); by=d.beta_outcome.values.astype(float); sy=d.se_outcome.values.astype(float); w=1/sy**2; den=np.sum(w*bx**2); b=np.sum(w*bx*by)/den; sf=np.sqrt(1/den); Q=np.sum(w*(by-b*bx)**2); df=len(d)-1; phi=max(1,Q/df) if random and df>0 else 1; se=sf*np.sqrt(phi); p=2*stats.norm.sf(abs(b/se)); return {'beta':b,'se':se,'p':p,'ci_low':b-1.96*se,'ci_high':b+1.96*se,'OR':math.exp(b),'OR_low':math.exp(b-1.96*se),'OR_high':math.exp(b+1.96*se),'Q':Q,'Q_df':df,'Q_p':stats.chi2.sf(Q,df) if df>0 else None,'phi':phi}

def egger(d):
 bx=d.beta.values.astype(float); by=d.beta_outcome.values.astype(float); sy=d.se_outcome.values.astype(float); fit=sm.WLS(by,sm.add_constant(bx),weights=1/sy**2).fit(); ci=fit.conf_int(); return {'intercept':fit.params[0],'intercept_se':fit.bse[0],'intercept_p':fit.pvalues[0],'beta':fit.params[1],'se':fit.bse[1],'p':fit.pvalues[1],'ci_low':ci[1,0],'ci_high':ci[1,1],'OR':math.exp(fit.params[1]),'OR_low':math.exp(ci[1,0]),'OR_high':math.exp(ci[1,1])}

def wm_point(d,by=None):
 bx=d.beta.values.astype(float); yy=d.beta_outcome.values.astype(float) if by is None else by; sy=d.se_outcome.values.astype(float); r=yy/bx; w=bx**2/sy**2; ix=np.argsort(r); r=r[ix]; w=w[ix]; c=np.cumsum(w)/w.sum(); j=np.searchsorted(c,.5); return float(r[j])

def wm(d,N=5000):
 b=wm_point(d); rng=np.random.default_rng(260908); by=d.beta_outcome.values.astype(float); sy=d.se_outcome.values.astype(float); sims=np.array([wm_point(d,rng.normal(by,sy)) for _ in range(N)]); se=sims.std(ddof=1); lo,hi=np.quantile(sims,[.025,.975]); return {'beta':b,'se':se,'p':2*stats.norm.sf(abs(b/se)),'ci_low':lo,'ci_high':hi,'OR':math.exp(b),'OR_low':math.exp(lo),'OR_high':math.exp(hi)}

def huber(d):
 bx=d.beta.values.astype(float); by=d.beta_outcome.values.astype(float); sy=d.se_outcome.values.astype(float); fit=sm.RLM(by/sy,(bx/sy)[:,None],M=sm.robust.norms.HuberT()).fit(); b=fit.params[0]; se=fit.bse[0]; return {'beta':b,'se':se,'p':2*stats.norm.sf(abs(b/se)),'ci_low':b-1.96*se,'ci_high':b+1.96*se,'OR':math.exp(b),'OR_low':math.exp(b-1.96*se),'OR_high':math.exp(b+1.96*se)}

def calc(d):
 loo=[]
 for i,r in d.iterrows():
  q=ivw(d.drop(i),True); loo.append({'left_out':r.rsid,'beta':q['beta'],'se':q['se'],'p':q['p']})
 l=pd.DataFrame(loo)
 return {'n':len(d),'mean_F':d.F.mean(),'min_F':d.F.min(),'max_F':d.F.max(),'ivw_mre':ivw(d,True),'ivw_fixed':ivw(d,False),'egger':egger(d),'weighted_median':wm(d),'robust_huber':huber(d),'loo_beta_min':l.beta.min(),'loo_beta_max':l.beta.max(),'loo_p_min':l.p.min(),'loo_p_max':l.p.max()},l

# female already harmonized from timed job
fd=pd.read_csv(B+'/female_harmonized_main.csv')
# male from hits
mi=pd.read_csv(B+'/main_instr/male_instruments.csv'); mo=pd.read_csv(B+'/male_outcome_hits.tsv',sep='\t'); md,mdrops=harmonize(mi,mo); md.to_csv(B+'/male_harmonized_main.csv',index=False)
fr,floo=calc(fd); mr,mloo=calc(md); floo.to_csv(B+'/female_loo_main.csv',index=False); mloo.to_csv(B+'/male_loo_main.csv',index=False)
f=fr['ivw_mre']; m=mr['ivw_mre']; diff=f['beta']-m['beta']; sed=math.sqrt(f['se']**2+m['se']**2); z=diff/sed; ip=2*stats.norm.sf(abs(z))
inter={'beta_difference_female_minus_male':diff,'se_difference':sed,'z':z,'p':ip,'ci_low':diff-1.96*sed,'ci_high':diff+1.96*sed}
wf=1/f['se']**2; wm_=1/m['se']**2; b=(wf*f['beta']+wm_*m['beta'])/(wf+wm_); se=math.sqrt(1/(wf+wm_)); meta={'beta':b,'se':se,'p':2*stats.norm.sf(abs(b/se)),'ci_low':b-1.96*se,'ci_high':b+1.96*se,'OR':math.exp(b),'OR_low':math.exp(b-1.96*se),'OR_high':math.exp(b+1.96*se)}
out={'female':fr,'male':mr,'male_drops':mdrops,'sex_interaction':inter,'fixed_effect_meta':meta}
open(B+'/main_mr_summary.json','w').write(json.dumps(out,indent=2,default=float))
print(json.dumps(out,indent=2,default=float))
print('\nFEMALE HARMONIZED RSIDs', ','.join(fd.rsid))
print('MALE HARMONIZED RSIDs', ','.join(md.rsid))
