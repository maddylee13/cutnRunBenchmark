##TSS_heatmap.sh is a script used to generate a heatmap over transcription start sites

#!/bin/bash
#PBS -l select=1:mem=128gb:ncpus=10
#PBS -l walltime=48:00:00
#PBS -N 241122_heatmap_CR_000

module load anaconda3/personal
source activate deeptools

DIR=/rds/general/user/mml120/ephemeral/p2/CR_000/alignment/bigwig
genes=/rds/general/user/mml120/ephemeral/p2/CR_000/genes
OUTDIR=/rds/general/user/mml120/ephemeral/p2/CR_000/heatmap

cores=8

computeMatrix scale-regions -S $DIR/CR_000_C2A_raw.bw \
                               $DIR/CR_000_C2B_raw.bw \
                               $DIR/CR_000_C5A_raw.bw \
                               $DIR/CR_000_C5B_raw.bw \
                               $DIR/CR_000_E2A_raw.bw \
                               $DIR/CR_000_E2B_raw.bw \
                               $DIR/CR_000_E5A_raw.bw \
                               $DIR/CR_000_E5B_raw.bw \
                               $DIR/CR_000_M2A_raw.bw \
                               $DIR/CR_000_M2B_raw.bw \
                               $DIR/CR_000_M5A_raw.bw \
                               $DIR/CR_000_M5B_raw.bw \
                              -R $genes/T2T-CHM13v2.0_genomic_modified.gtf \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 5000 \
                              --afterRegionStartLength 3000 \
                              --skipZeros -o $OUTDIR/T2T_gene/CR_000_matrix_gene.mat.gz -p $cores

plotHeatmap -m $OUTDIR/T2T_gene/CR_000_matrix_gene.mat.gz -out $OUTDIR/T2T_gene/CT_000_Histone_gene.png --sortUsing sum
