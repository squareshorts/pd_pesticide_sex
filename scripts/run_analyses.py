from pathlib import Path
import concurrent.futures,subprocess
B=Path(__file__).resolve().parents[1]
jobs=[f'{c}_{s}' for c in ['22614_1','22614_0','22614_2','22610_2'] for s in ['female','male']]
def run(key):
 if (B/'results'/f'{key}_analysis.json').exists(): return 0
 with (B/'logs'/f'{key}_analysis.log').open('w') as f:
  r=subprocess.run(['C:/Program Files/R/R-4.5.2/bin/Rscript.exe',str(B/'scripts/analyze.R'),key],cwd=B,stdout=f,stderr=subprocess.STDOUT)
 print(key,r.returncode,flush=True)
 return r.returncode
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 status=list(pool.map(run,jobs))
assert not any(status),status
