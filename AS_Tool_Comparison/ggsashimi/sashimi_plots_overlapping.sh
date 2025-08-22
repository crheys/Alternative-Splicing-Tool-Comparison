#!/bin/sh

#event - 22:42601969-42602056 
#POLDIP3 (SE)
/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/ggsashimi.py \
  -b input_bams_fixed.tsv \
  -c 22:42599740-42603000 \
  -g /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/reference/GRCh38.p14/Homo_sapiens.GRCh38.113.gtf \
  -M 10 -C 3 -O 3 \
  --shrink --alpha 0.25 \
  --base-size=12 --ann-height=8 --height=6 --width=18 \
  --fix-y-scale \
  -P palette.txt \
   -A mean \
  -o POLDIP3_FTD_se

#event - 19:17641556:17642414
#UNC13A (cryptic exon)
/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/ggsashimi.py \
  -b input_bams_fixed.tsv \
  -c 19:17641400-17642950 \
  -g /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/reference/GRCh38.p14/Homo_sapiens.GRCh38.113.gtf \
  -M 10 -C 3 -O 3 \
  --shrink --alpha 0.25 \
  --base-size=12 --ann-height=8 --height=6 --width=18 \
  --fix-y-scale \
  -P palette.txt \
   -A mean \
  -o UNC13A_FTD

#event - 12:88086498-88086642
#CEP290 (cryptic exon)
/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/ggsashimi.py \
  -b input_bams_fixed.tsv \
  -c 12:88086410-88087900 \
  -g /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/reference/GRCh38.p14/Homo_sapiens.GRCh38.113.gtf \
  -M 10 -C 3 -O 3 \
  --shrink --alpha 0.25 \
  --base-size=12 --ann-height=8 --height=6 --width=18 \
  --fix-y-scale \
  -P palette.txt \
   -A mean \
  -o CEP290_leafcutter
