import pathlib, gzip, json, tarfile, zipfile, subprocess, hashlib
import pandas as pd
import numpy as np
BASE=pathlib.Path(__file__).resolve().parents[1]; D=BASE/'data'; O=BASE/'results'
codes=['22614_0','22614_1','22614_2','22610_2']; sexes=['female','male']
def extract_candidates():
 for code in codes:
  for sex in sexes:
   key=f'{code}_{sex}'; dest=D/f'{key}_candidates_raw.csv'
   if dest.exists(): continue
   path=D/f'{code}.gwas.imputed_v3.{sex}.tsv.bgz'
   if not path.exists() or path.stat().st_size<500000000:
    print('PENDING',key,flush=True); continue
   rows=[]; nall=0; nvalid=0; nsig=0; nconf=0; nmaf=0; nf=0
   for x in pd.read_csv(path,sep='\t',compression='gzip',chunksize=250000):
    nall+=len(x)
    x=x[np.isfinite(x.beta)&np.isfinite(x.se)&(x.se>0)&np.isfinite(x.pval)]; nvalid+=len(x)
    x=x[x.pval<1e-5]; nsig+=len(x)
    x=x[x.low_confidence_variant.astype(str).str.lower()=='false']; nconf+=len(x)
    x=x[x.minor_AF>=.01]; nmaf+=len(x)
    x['F']=(x.beta/x.se)**2; x=x[x.F>10]; nf+=len(x)
    rows.append(x)
   x=pd.concat(rows,ignore_index=True)
   v=x.variant.str.split(':',expand=True); x['chrom']=v[0]; x['pos']=v[1].astype(int); x['other_allele']=v[2]; x['effect_allele']=v[3]
   x['eaf']=x.AC/(2*x.n_complete_samples); x['n']=x.n_complete_samples
   x.to_csv(dest,index=False)
   counts=dict(total=nall,valid_numeric=nvalid,p_lt_1e5=nsig,high_confidence=nconf,MAF_ge_001=nmaf,F_gt_10=nf)
   (D/f'{key}_stages.json').write_text(json.dumps(counts,indent=2)); print(key,counts,flush=True)
def reference():
 if not (D/'EUR.bim').exists():
  with tarfile.open(D/'1kg.v3.tgz','r:gz') as t:
   for m in t:
    if pathlib.PurePosixPath(m.name).name in ['EUR.bed','EUR.bim','EUR.fam']:
     out=D/pathlib.PurePosixPath(m.name).name
     with t.extractfile(m) as src, out.open('wb') as dst:
      import shutil; shutil.copyfileobj(src,dst)
 if not (D/'outcome_files.json').exists():
  outer=D/'Blauwendraat_Parkinsons_sexdiff_2021.zip'; inner=D/'Blauwendraat_male_female_GWAS.zip'
  if not inner.exists():
   with zipfile.ZipFile(outer) as z: z.extract('Blauwendraat_male_female_GWAS.zip',D)
  with zipfile.ZipFile(inner) as z:
   names=[n for n in z.namelist() if 'NO_UKB_AT_ALL_no_multi_allelics_RSID.txt.gz' in n and not n.startswith('__MACOSX')]
   print('Outcome entries',names,flush=True)
   for n in names:
    with z.open(n) as src,(D/pathlib.PurePosixPath(n).name).open('wb') as dst:
     import shutil; shutil.copyfileobj(src,dst)
  (D/'outcome_files.json').write_text(json.dumps(names,indent=2))
def clump():
 # Map positions AND unordered SNP alleles against European reference; no position-only mapping.
 wanted=set()
 for path in D.glob('*_candidates_raw.csv'):
  x=pd.read_csv(path); wanted.update(zip(x.chrom.astype(str),x.pos))
 mapping={}
 with (D/'EUR.bim').open() as f:
  for line in f:
   ch,rs,cm,pos,a,b=line.split(); pos=int(pos)
   if (ch,pos) in wanted:
    mapping.setdefault((ch,pos,tuple(sorted([a,b]))),[]).append(rs)
 for code in codes:
  for sex in sexes:
   key=f'{code}_{sex}'; x=pd.read_csv(D/f'{key}_candidates_raw.csv'); counts=json.loads((D/f'{key}_stages.json').read_text())
   def maprs(r):
    m=mapping.get((str(r.chrom),r.pos,tuple(sorted([r.effect_allele,r.other_allele]))),[])
    return m[0] if len(m)==1 else None
   x['rsid']=x.apply(maprs,axis=1); x.to_csv(O/f'{key}_candidates.csv',index=False)
   x=x.dropna(subset=['rsid']).sort_values(['pval','rsid']).drop_duplicates('rsid'); counts['reference_matched_unique']=len(x)
   inp=O/f'{key}_clump_input.tsv'; x[['rsid','pval']].rename(columns={'rsid':'SNP','pval':'P'}).to_csv(inp,sep='\t',index=False)
   prefix=O/f'{key}_ld'
   cmd=[str(BASE/'software/plink/plink.exe'),'--bfile',str(D/'EUR'),'--clump',str(inp),'--clump-p1','0.00001','--clump-p2','0.00001','--clump-r2','0.001','--clump-kb','10000','--threads','2','--memory','2048','--out',str(prefix)]
   subprocess.run(cmd,check=True,stdout=(BASE/'logs'/f'{key}_plink_stdout.txt').open('w'))
   leads=pd.read_csv(str(prefix)+'.clumped',sep=r'\s+'); x=x[x.rsid.isin(leads.SNP)]; counts['clumped']=len(x)
   x.to_csv(O/f'{key}_instruments.csv',index=False); (O/f'{key}_stages.json').write_text(json.dumps(counts,indent=2)); print(key,'clumped',len(x),flush=True)
def harmonize():
 comp=str.maketrans('ACGT','TGCA')
 for sex in sexes:
  wanted=set()
  for code in codes: wanted.update(pd.read_csv(O/f'{code}_{sex}_instruments.csv').rsid)
  path=D/f'{sex.upper()}_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt.gz'
  found=[]
  for x in pd.read_csv(path,sep='\t',compression='gzip',chunksize=250000):
   found.append(x[x.ID.isin(wanted)])
  out=pd.concat(found,ignore_index=True); out.to_csv(O/f'{sex}_outcome_hits.csv',index=False)
  assert not out.ID.duplicated().any()
  om=out.set_index('ID').to_dict('index')
  for code in codes:
   key=f'{code}_{sex}'; inst=pd.read_csv(O/f'{key}_instruments.csv'); counts=json.loads((O/f'{key}_stages.json').read_text()); rows=[]; audit=[]
   for _,r in inst.iterrows():
    rs=r.rsid; a=r.effect_allele.upper(); b=r.other_allele.upper(); pal={a,b} in [{'A','T'},{'C','G'}]
    rec={'rsid':rs,'palindromic':pal,'status':'retained'}
    if rs not in om: rec['status']='absent_outcome'; audit.append(rec); continue
    y=om[rs]; c=str(y['Allele1']).upper(); e=str(y['Allele2']).upper(); fy=float(y['Freq1']); by=float(y['Effect']); sy=float(y['StdErr'])
    flip=None; strand=False
    for ac,bc,st in [(a,b,False),(a.translate(comp),b.translate(comp),True)]:
     if (ac,bc)==(c,e): flip=False; strand=st; break
     if (ac,bc)==(e,c): flip=True; strand=st; break
    if flip is None: rec['status']='allele_mismatch'
    elif not np.isfinite([fy,by,sy]).all() or sy<=0 or fy<=0 or fy>=1: rec['status']='invalid_outcome'
    elif pal:
     if (.42<=r.eaf<=.58) or (.42<=fy<=.58): rec['status']='palindromic_ambiguous_frequency'
     else:
      af=1-fy if flip else fy
      if (r.eaf<.5)!=(af<.5): flip=not flip
      if abs(r.eaf-(1-fy if flip else fy))>.20: rec['status']='palindromic_frequency_mismatch'
    if rec['status']=='retained':
     af=1-fy if flip else fy
     rec.update(flip_outcome=flip,strand_complement=strand,frequency_difference=abs(r.eaf-af))
     rows.append({**r.to_dict(),'beta_outcome':-by if flip else by,'se_outcome':sy,'p_outcome':float(y['P-value']),'eaf_outcome_aligned':af,'outcome_effect_allele_raw':c,'outcome_other_allele_raw':e,'flip_outcome':flip,'strand_complement':strand,'palindromic':pal})
    audit.append(rec)
   dat=pd.DataFrame(rows); dat.to_csv(O/f'{key}_harmonized.csv',index=False); pd.DataFrame(audit).to_csv(O/f'{key}_harmonization_audit.csv',index=False)
   counts['outcome_found']=sum(r in om for r in inst.rsid); counts['harmonized']=len(dat); counts['harmonization_drops']=pd.DataFrame(audit).status.value_counts().to_dict()
   (O/f'{key}_stages.json').write_text(json.dumps(counts,indent=2)); print(key,counts,flush=True)
if __name__=='__main__':
 import sys
 if sys.argv[-1]=='candidates': extract_candidates()
 else: reference(); clump(); harmonize()
