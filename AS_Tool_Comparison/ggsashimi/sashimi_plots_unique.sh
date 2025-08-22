#!/bin/sh

#event - 14:60812127-60818722
#MNAT1 (cryptic exon)
/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/ggsashimi.py \
  -b input_bams_fixed.tsv \
  -c 14:60812000-60818825 \
  -g /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/reference/GRCh38.p14/Homo_sapiens.GRCh38.113.gtf \
  -M 10 -C 3 -O 3 \
  --shrink --alpha 0.25 \
  --base-size=12 --ann-height=8 --height=6 --width=18 \
  --fix-y-scale \
  -P palette.txt \
   -A mean \
  -o MNAT1_leafcutter

#event - 22:31845018-31845237 
#DEPDC5 SE 
/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/ggsashimi.py \
  -b input_bams_fixed.tsv \
  -c 22:31843650-31846950 \
  -g /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/reference/GRCh38.p14/Homo_sapiens.GRCh38.113.gtf \
  -M 1 -C 3 -O 3 \
  --shrink --alpha 0.25 \
  --base-size=12 --ann-height=8 --height=6 --width=18 \
  --fix-y-scale \
  -P palette.txt \
   -A mean \
  -o DEPDC5_DEXSeq

#event - X:102713091-102713217
#GPRASP2 SE 
/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/ggsashimi.py \
  -b input_bams_fixed.tsv \
  -c X:102712500-102713830 \
  -g /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/reference/GRCh38.p14/Homo_sapiens.GRCh38.113.gtf \
  -M 1 -C 3 -O 3 \
  --shrink --alpha 0.25 \
  --base-size=12 --ann-height=8 --height=6 --width=18 \
  --fix-y-scale \
  -P palette.txt \
   -A mean \
  -o GPRASP2_SUPPA2

#event - 8:81717249-81717288
#ZFAND1 SE
/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/ggsashimi.py \
  -b input_bams_fixed.tsv \
  -c 8:81715000-81718220 \
  -g /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/reference/GRCh38.p14/Homo_sapiens.GRCh38.113.gtf \
  -M 10 -C 3 -O 3 \
 --shrink --alpha 0.25 \
  --base-size=12 --ann-height=8 --height=6 --width=18 \
  --fix-y-scale \
  -P palette.txt \
   -A mean \
  -o ZFAND1_MAJIQ

#event - 16:8929368-8929543
#USP7 SE 
/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/ggsashimi.py \
  -b input_bams_fixed.tsv \
  -c 16:8923214-8930397 \
  -g /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/reference/GRCh38.p14/Homo_sapiens.GRCh38.113.gtf \
  -M 1 -C 3 -O 3 \
  --shrink --alpha 0.25 \
  --base-size=12 --ann-height=8 --height=6 --width=18 \
  --fix-y-scale \
  -P palette.txt \
   -A mean \
  -o USP7_rMATS

#event - 11:66423714-66423834
#NPAS4 RI 
/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/ggsashimi.py \
  -b input_bams_fixed.tsv \
  -c 11:66423600-66424100 \
  -g /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/reference/GRCh38.p14/Homo_sapiens.GRCh38.113.gtf \
  -M 3 -C 3 -O 3 \
--alpha 0.25 \
  --base-size=12 --ann-height=8 --height=6 --width=18 \
  --fix-y-scale \
  -P palette.txt \
   -A mean \
  -o NPAS4_whippet
  
  
