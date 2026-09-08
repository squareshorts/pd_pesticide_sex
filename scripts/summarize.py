import pathlib,json,math,hashlib,platform
import numpy as np,pandas as pd,scipy
from scipy import stats
B=pathlib.Path(__file__).resolve().parents[1]; O=B/'results'
codes=['22614_1','22614_0','22614_2','22610_2']; sexes=['female','male']; zcrit=stats.norm.ppf(.975)
def normresult(b,se):
 lo=b-zcrit*se; hi=b+zcrit*se
 return dict(beta=b,se=se,p=2*stats.norm.sf(abs(b/se)),ci_low=lo,ci_high=hi,OR=np.exp(b),OR_low=np.exp(lo),OR_high=np.exp(hi))
def bh(p):
 p=np.asarray(p); ix=np.argsort(p); q=np.empty(len(p)); q[ix]=np.minimum(1,np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]); return q
def median_point(x,y,w):
 ix=np.argsort(y/x); r=(y/x)[ix]; w=w[ix]; c=(np.cumsum(w)-w/2)/sum(w)
 return np.interp(.5,c,r)
all={}; core=[]; methods=[]; sens=[]; stages=[]; verified=[]
for code in codes:
 for sex in sexes:
  key=f'{code}_{sex}'; a=json.loads((O/f'{key}_analysis.json').read_text()); d=pd.read_csv(O/f'{key}_harmonized.csv'); all[key]=a
  x=d.beta.values; y=d.beta_outcome.values; sx=d.se.values; sy=d.se_outcome.values
  den=np.sum(x*x/(sy*sy)); b=np.sum(x*y/(sy*sy))/den; q=np.sum((y-b*x)**2/(sy*sy)); se=np.sqrt(max(1,q/(len(x)-1))/den)
  for k,v in normresult(b,se).items(): assert np.isclose(v,a['ivw'][k],rtol=1e-10,atol=1e-12),(key,k,v,a['ivw'][k])
  assert np.isclose(q,a['ivw']['Q'])
  for method,w in [('weighted_median',x*x/(sy*sy)),('weighted_median_second_order',1/(sy**2/x**2+y**2*sx**2/x**4))]:
   assert np.isclose(median_point(x,y,w),a[method]['beta'],atol=1e-10)
  X=np.column_stack([np.ones(len(x)),abs(x)]); Y=y*np.sign(x); w=1/sy**2
  inv=np.linalg.inv(X.T@(w[:,None]*X)); coef=inv@(X.T@(w*Y)); qe=np.sum(w*(Y-X@coef)**2); ses=np.sqrt(np.diag(inv)*max(1,qe/(len(x)-2)))
  for k,v in [('beta',coef[1]),('se',ses[1]),('intercept',coef[0]),('intercept_se',ses[0]),('p',2*stats.t.sf(abs(coef[1]/ses[1]),len(x)-2)),('intercept_p',2*stats.t.sf(abs(coef[0]/ses[0]),len(x)-2))]: assert np.isclose(v,a['egger'][k],rtol=1e-9,atol=1e-11),(key,k)
  od=pd.read_csv(O/f'{key}_egger_oriented.csv'); flip=x<0
  assert np.allclose(od.beta,abs(x)) and np.allclose(od.beta_outcome,y*np.sign(x))
  assert np.allclose(od.se,sx) and np.allclose(od.se_outcome,sy)
  assert np.all(od.effect_allele==np.where(flip,d.other_allele,d.effect_allele))
  assert np.allclose(od.eaf,np.where(flip,1-d.eaf,d.eaf))
  qgx=np.sum((abs(x)-np.average(abs(x),weights=1/sx**2))**2/sx**2)
  assert np.isclose(max(0,1-(len(x)-1)/qgx),a['I2GX'])
  for loss in ['huber','l2']:
   rr=a['raps_'+loss]
   if 'beta.hat' in rr: a['raps_'+loss+'_result']=normresult(rr['beta.hat'],rr['beta.se'])
  e=a['egger']; e['OR']=np.exp(e['beta']); e['OR_low']=np.exp(e['ci_low']); e['OR_high']=np.exp(e['ci_high'])
  base=dict(exposure=code,sex=sex,n=a['n'])
  core.append({**base,'mean_F':a['F']['mean'],**a['ivw'],'WM_beta':a['weighted_median']['beta'],'WM_se':a['weighted_median']['se'],'WM_p':a['weighted_median']['p'],'egger_beta':e['beta'],'egger_se':e['se'],'egger_p':e['p'],'egger_intercept':e['intercept'],'egger_intercept_se':e['intercept_se'],'egger_intercept_p':e['intercept_p'],'I2GX':a['I2GX'],**{'RAPS_'+k:v for k,v in a.get('raps_huber_result',{}).items()}})
  for method in ['ivw','weighted_median','weighted_median_second_order','egger','raps_huber_result','raps_l2_result']:
   if method in a: methods.append({**base,'method':method,**a[method]})
  mo=a['weighted_mode']; mo=mo[0] if isinstance(mo,list) else mo
  if 'b' in mo: methods.append({**base,'method':'weighted_mode','beta':mo['b'],'se':mo['se'],'p':mo['pval'],'ci_low':mo.get('ci_low'),'ci_high':mo.get('ci_upp')})
  pr=a.get('presso',{}); globalp=pr.get('MR-PRESSO results',{}).get('Global Test',{}).get('Pvalue')
  raw=pr.get('Main MR results',[])
  if isinstance(raw,list):
   for r in raw:
    if r.get('Causal Estimate') is not None: methods.append({**base,'method':'MR_PRESSO_'+r['MR Analysis'],'beta':r['Causal Estimate'],'se':r['Sd'],'p':r['P-value']})
  sens.append({**base,**{'F_'+k:v for k,v in a['F'].items()},'Q':a['ivw']['Q'],'Q_df':a['ivw']['Q_df'],'Q_p':a['ivw']['Q_p'],'I2GX':a['I2GX'],'QGX':a['QGX'],'egger_intercept':e['intercept'],'egger_intercept_se':e['intercept_se'],'egger_intercept_p':e['intercept_p'],**{'loo_'+k:v for k,v in a['loo'].items()},'radial_n_outliers':a['radial']['n_outliers'],'PRESSO_global_p':globalp,'AR_p_null':a['anderson_rubin']['p_null'],'AR_95_confidence_set':json.dumps(a['anderson_rubin']['confidence_set']),'palindromic_retained':int(d.palindromic.sum()),'n_egger_flipped':int(flip.sum()),'steiger_status':a['steiger']['status']})
  st=json.loads((O/f'{key}_stages.json').read_text()); st['harmonization_drops']=json.dumps(st['harmonization_drops']); stages.append({**base,**st})
  verified.append({'analysis':key,'independent_IVW':'pass','independent_Egger':'pass','weighted_median_points':'pass','I2GX':'pass','allele_reorientation':'pass'})
pd.DataFrame(verified).to_csv(O/'numerical_validation.csv',index=False)
inter=[]; pooled=[]
for code in codes:
 f=all[code+'_female']['ivw']; m=all[code+'_male']['ivw']; d=f['beta']-m['beta']; se=math.hypot(f['se'],m['se']); z=d/se
 inter.append(dict(exposure=code,test_role='primary' if code=='22614_1' else 'secondary',beta_difference_female_minus_male=d,se_difference=se,Z=z,p=2*stats.norm.sf(abs(z)),ci_low=d-zcrit*se,ci_high=d+zcrit*se))
 w=np.array([1/f['se']**2,1/m['se']**2]); b=np.average([f['beta'],m['beta']],weights=w); se=1/np.sqrt(sum(w))
 pooled.append(dict(exposure=code,role='secondary_descriptive',**normresult(b,se)))
inter=pd.DataFrame(inter); inter['FDR_BH_secondary_3']=np.nan; ix=inter.exposure!='22614_1'; inter.loc[ix,'FDR_BH_secondary_3']=bh(inter.loc[ix,'p'])
core=pd.DataFrame(core); core['FDR_BH_secondary_IVW_6']=np.nan; ix=core.exposure!='22614_1'; core.loc[ix,'FDR_BH_secondary_IVW_6']=bh(core.loc[ix,'p'])
# Sensitivity to defining one wider family: all six secondary IVWs plus three interactions.
q9=bh(list(core.loc[ix,'p'])+list(inter.loc[inter.exposure!='22614_1','p'])); core['FDR_BH_secondary_all_9']=np.nan; core.loc[ix,'FDR_BH_secondary_all_9']=q9[:6]; inter['FDR_BH_secondary_all_9']=np.nan; inter.loc[inter.exposure!='22614_1','FDR_BH_secondary_all_9']=q9[6:]
tables={'primary_results':core[core.exposure=='22614_1'],'all_IVW_results':core,'sex_interactions':inter,'pooled_descriptive':pd.DataFrame(pooled),'all_estimators':pd.DataFrame(methods),'sensitivity_summary':pd.DataFrame(sens),'filtering_stages':pd.DataFrame(stages)}
secondary=core[core.exposure!='22614_1'].merge(inter,on='exposure',suffixes=('_IVW','_interaction')); tables['secondary_panel']=secondary
for name,df in tables.items(): df.to_csv(O/f'{name}.csv',index=False)
(O/'all_analyses.json').write_text(json.dumps(all,indent=2,allow_nan=False))
# A separate, recomputed audit isolates the earlier allele-coding and weighted-quantile problems.
saved=json.loads((O/'saved_primary_correction_audit.json').read_text()); rows=[]
for sex in sexes:
 d=pd.read_csv(B/f'{sex}_harmonized_main.csv'); r=d.beta_outcome.values/d.beta.values; w=d.beta.values**2/d.se_outcome.values**2; ix=np.argsort(r); oldwm=r[ix][np.searchsorted(np.cumsum(w[ix])/sum(w),.5)]
 a=saved[sex]
 rows.append(dict(sex=sex,old_egger_beta_recomputed=a['previous_algorithm_recomputed']['beta'],old_egger_se_recomputed=a['previous_algorithm_recomputed']['se'],old_egger_p_recomputed=a['previous_algorithm_recomputed']['p'],old_intercept_p_recomputed=a['previous_algorithm_recomputed']['intercept_p'],orientation_only_beta=a['orientation_only_recomputed']['beta'],orientation_only_se=a['orientation_only_recomputed']['se'],orientation_only_p=a['orientation_only_recomputed']['p'],corrected_egger_beta=a['corrected_egger']['beta'],corrected_egger_se=a['corrected_egger']['se'],corrected_egger_p=a['corrected_egger']['p'],corrected_intercept_p=a['corrected_egger']['intercept_p'],old_weighted_median_point_recomputed=oldwm,corrected_weighted_median_beta=a['analysis']['weighted_median']['beta'],corrected_weighted_median_p=a['analysis']['weighted_median']['p'],saved_ivw_beta_recomputed=a['analysis']['ivw']['beta'],new_non_v2_ivw_beta=all['22614_1_'+sex]['ivw']['beta']))
pd.DataFrame(rows).to_csv(O/'earlier_results_changes.csv',index=False)
print('PRIMARY\n',tables['primary_results'].to_string(index=False)); print('INTERACTIONS\n',inter.to_string(index=False)); print('SENSITIVITIES\n',tables['sensitivity_summary'].to_string(index=False)); print('POOLED\n',tables['pooled_descriptive'].to_string(index=False))
(O/'python_versions.json').write_text(json.dumps(dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__),indent=2))

