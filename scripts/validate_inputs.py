import pathlib,subprocess,hashlib,gzip,json,concurrent.futures
import numpy as np,pandas as pd
B=pathlib.Path(__file__).resolve().parents[1]; D=B/'data'; O=B/'results'
checks=[]
for path in sorted(O.glob('226*_harmonized.csv')):
 d=pd.read_csv(path); key=path.stem.replace('_harmonized',''); ref=pd.read_csv(O/f'{key}_instruments.csv'); out=pd.read_csv(O/f'{key.split("_")[-1]}_outcome_hits.csv').set_index('ID')
 assert d.rsid.is_unique and d.rsid.isin(ref.rsid).all()
 assert (d.pval<1e-5).all() and (d.F>10).all() and (d.minor_AF>=.01).all()
 assert np.allclose(d.F,(d.beta/d.se)**2)
 assert np.allclose(np.minimum(d.eaf,1-d.eaf),d.minor_AF,atol=1e-5)
 for _,r in d.iterrows():
  y=out.loc[r.rsid]; flip=str(r.flip_outcome).lower()=='true'
  assert np.isclose(r.beta_outcome,(-1 if flip else 1)*y.Effect)
  assert r.se_outcome==y.StdErr
  assert np.isclose(r.eaf_outcome_aligned,1-y.Freq1 if flip else y.Freq1)
 ids=O/f'{key}_retained_ids.txt'; ids.write_text('\n'.join(d.rsid)+'\n')
 prefix=O/f'{key}_retained_ld_check'
 cmd=[str(B/'software/plink/plink.exe'),'--bfile',str(D/'EUR'),'--extract',str(ids),'--r2','dprime','--ld-window','999999','--ld-window-kb','10000','--ld-window-r2','0.001','--memory','2048','--threads','2','--out',str(prefix)]
 subprocess.run(cmd,check=True,stdout=(B/'logs'/f'{key}_ldcheck_stdout.txt').open('w'))
 ld=pd.read_csv(str(prefix)+'.ld',sep=r'\s+')
 assert len(ld)==0,(key,ld)
 checks.append({'analysis':key,'n':len(d),'selection_harmonization_checks':'pass','retained_pairs_r2_ge_001_within_10Mb':len(ld)})
 print('PASS',key,flush=True)
pd.DataFrame(checks).to_csv(O/'data_validation.csv',index=False)
def outcome_hash(sex):
 name=f'{sex}_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt.gz'; path=D/name
 expected={'FEMALE':'4c3123aa1622fac6653efe41662faeb1','MALE':'2a400dbbda127da0635f3ef72f36421b'}[sex]
 with gzip.open(path,'rb') as f: md5=hashlib.file_digest(f,'md5').hexdigest()
 assert md5==expected,(name,md5,expected)
 return dict(file=name,uncompressed_md5=md5,expected_uncompressed_md5=expected,sha256=hashlib.file_digest(path.open('rb'),'sha256').hexdigest())
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool: hashes=list(pool.map(outcome_hash,['FEMALE','MALE']))
(O/'outcome_checksums.json').write_text(json.dumps(hashes,indent=2))
print('Outcome uncompressed MD5 checks passed',flush=True)

