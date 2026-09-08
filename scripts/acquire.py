import concurrent.futures, subprocess, pathlib, csv, json, hashlib
BASE=pathlib.Path(__file__).resolve().parents[1]
D=BASE/'data'
jobs=[]
for code in ['22614_0','22614_1','22614_2','22610_2']:
 for sex in ['female','male']:
  name=f'{code}.gwas.imputed_v3.{sex}.tsv.bgz'
  jobs.append((name,'https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/'+name))
jobs += [('Blauwendraat_Parkinsons_sexdiff_2021.zip','https://personal.broadinstitute.org/ryank/Blauwendraat_Parkinsons_sexdiff_2021.zip'),('1kg.v3.tgz','http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz')]
def get(job):
 name,url=job; path=D/name
 if path.exists() and path.stat().st_size>500000000 and name.endswith('.bgz'):
  return {'file':name,'url':url,'returncode':0,'bytes':path.stat().st_size,'sha256':hashlib.file_digest(path.open('rb'),'sha256').hexdigest()}
 print('START',name,flush=True)
 r=subprocess.run(['C:/Windows/System32/curl.exe','-L','--fail','--silent','--show-error','--retry','3','-C','-','--connect-timeout','30','--speed-limit','1024','--speed-time','60',url,'-o',str(path)],capture_output=True,text=True)
 result={'file':name,'url':url,'returncode':r.returncode,'stderr':r.stderr,'bytes':path.stat().st_size if path.exists() else 0}
 if r.returncode==0:
  result['sha256']=hashlib.file_digest(path.open('rb'),'sha256').hexdigest()
 print('DONE',result,flush=True)
 return result
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
 results=list(pool.map(get,jobs))
(D/'download_manifest.json').write_text(json.dumps(results,indent=2))
