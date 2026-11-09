"""EUR LD proxies (r2 >= 0.8, +/-1 Mb) for the 37 primary instruments, from data/EUR.*; one window at a time."""
import numpy as np, pandas as pd, ld_tools as L
p=pd.read_csv('primary_panel_37.csv'); out=[]
for _,r in p.iterrows():
    win=L.bim_lookup([],[(str(r.chrom),int(r.pos)-1_000_000,int(r.pos)+1_000_000)])
    G=L.genotypes([w[0] for w in win]); t=[w for w in win if w[2]==r.rsid][0]; gt=G[t[0]]
    for i,c,rs,pos,a1,a2 in win:
        if rs=='.' : continue
        v=L.r2(gt,G[i])
        if v>=0.8: out.append((r.rsid,r.origin,rs,c,pos,round(v,4)))
    del G
df=pd.DataFrame(out,columns=['instrument','origin','lookup_rsid','chr','pos_GRCh37','r2']); df.to_csv('eur_proxies_r2_0.8.csv',index=False)
print(df.groupby('instrument').size().describe()); print('unique lookup rsIDs',df.lookup_rsid.nunique()); print(df[df.instrument=='rs62092069'].shape)
