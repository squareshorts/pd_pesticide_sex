#!/bin/bash
# Extract the 37 primary variants from both sex-specific Neale 22614_1 files and both IPDGC outcome files.
cd "$(dirname "$0")"
D=../data
for s in female male; do
  zcat $D/22614_1.gwas.imputed_v3.$s.tsv.bgz | awk -F'\t' 'NR==FNR{w[$1]=1;next} FNR==1||($1 in w)' variants.txt - > neale_22614_1_${s}_37.tsv &
done
for S in FEMALE MALE; do
  zcat $D/${S}_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt.gz | awk -F'\t' 'NR==FNR{w[$1]=1;next} FNR==1||($16 in w)' rsids.txt - > ipdgc_${S}_37.tsv &
done
wait; echo DONE > extract.done
