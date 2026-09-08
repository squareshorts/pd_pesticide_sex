import pathlib,sys,concurrent.futures,hashlib,json,time
B=pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0,str(B/'software/pylib'))
import requests
D=B/'data'; chunks=D/'download_chunks'; chunks.mkdir(exist_ok=True)
jobs=[('22614_0.gwas.imputed_v3.male.tsv.bgz','https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/22614_0.gwas.imputed_v3.male.tsv.bgz',589991522),('1kg.v3.tgz','http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz',1567805147)]
def piece(args):
 name,url,start,end=args; p=chunks/f'{name}.{start}'
 if p.exists() and p.stat().st_size==end-start+1: return p
 for attempt in range(4):
  try:
   r=requests.get(url,headers={'Range':f'bytes={start}-{end}'},timeout=(20,60))
   r.raise_for_status(); assert r.status_code==206,(r.status_code,r.headers)
   assert r.headers['Content-Range'].startswith(f'bytes {start}-{end}/')
   assert len(r.content)==end-start+1
   p.write_bytes(r.content); print(name,start,end,'done',flush=True); return p
  except Exception as e:
   print(name,start,'retry',attempt,repr(e),flush=True)
 raise RuntimeError(args)
tasks=[]; offsets={}
for name,url,size in jobs:
 start=(D/name).stat().st_size if (D/name).exists() else 0; offsets[name]=start
 for a in range(start,size,16*1024*1024): tasks.append((name,url,a,min(size-1,a+16*1024*1024-1)))
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool: list(pool.map(piece,tasks))
for name,url,size in jobs:
 with (D/name).open('ab') as out:
  for a in range(offsets[name],size,16*1024*1024): out.write((chunks/f'{name}.{a}').read_bytes())
 assert (D/name).stat().st_size==size
 print('COMPLETE',name,hashlib.file_digest((D/name).open('rb'),'sha256').hexdigest(),flush=True)
