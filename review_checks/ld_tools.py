"""LD from the 1000 Genomes EUR PLINK binary set used in the original pipeline (data/EUR.*), without PLINK.
Genotypes are read directly from the SNP-major .bed file; r^2 is the squared Pearson correlation of
allele dosages over individuals with both genotypes observed (as in PLINK --r2)."""
import subprocess, numpy as np
from pathlib import Path
D=Path(__file__).resolve().parents[1]/'data'
N=sum(1 for _ in open(D/'EUR.fam')); BPS=(N+3)//4
_LUT=np.array([2,np.nan,1,0],dtype=float)  # PLINK codes 00,01,10,11 -> A1 dosage

def bim_lookup(rsids=(), windows=()):
    """Return rows (idx, chr, rsid, pos, a1, a2) for listed rsIDs and/or (chr, lo, hi) windows."""
    rs=' '.join(rsids); win=';'.join(f'{c},{lo},{hi}' for c,lo,hi in windows)
    prog=r'''BEGIN{n=split(RS_,a," ");for(i=1;i<=n;i++)w[a[i]]=1;m=split(WIN,b,";");for(j=1;j<=m;j++){split(b[j],c,",");C[j]=c[1];L[j]=c[2]+0;H[j]=c[3]+0}}
    {k=($2 in w);if(!k)for(j=1;j<=m;j++)if($1==C[j]&&$4>=L[j]&&$4<=H[j]){k=1;break} if(k)print NR-1"\t"$1"\t"$2"\t"$4"\t"$5"\t"$6}'''
    out=subprocess.run(['awk','-v',f'RS_={rs}','-v',f'WIN={win}',prog,str(D/'EUR.bim')],capture_output=True,text=True,check=True).stdout
    rows=[l.split('\t') for l in out.strip().split('\n') if l]
    return [(int(i),c,r,int(p),a,b) for i,c,r,p,a,b in rows]

def genotypes(indices):
    idx=np.asarray(sorted(set(indices)),dtype=np.int64); G=np.empty((len(idx),N))
    with open(D/'EUR.bed','rb') as f:
        assert f.read(3)==b'\x6c\x1b\x01'
        for k,i in enumerate(idx):
            f.seek(3+i*BPS); by=np.frombuffer(f.read(BPS),dtype=np.uint8)
            codes=np.stack([(by>>s)&3 for s in (0,2,4,6)],axis=1).ravel()[:N]
            G[k]=_LUT[codes]
    return dict(zip(idx.tolist(),G))

def r2(g1,g2):
    ok=~(np.isnan(g1)|np.isnan(g2))
    if ok.sum()<10 or g1[ok].std()==0 or g2[ok].std()==0: return np.nan
    return float(np.corrcoef(g1[ok],g2[ok])[0,1]**2)

def signed_r(g1,g2):
    ok=~(np.isnan(g1)|np.isnan(g2)); return float(np.corrcoef(g1[ok],g2[ok])[0,1])
