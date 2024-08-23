##macs2_Ig_heatmap.sh is a script used to generate heatmaps for MACS2 CUT&Run peaks with IgG as a control

#!/bin/bash
#PBS -l select=1:mem=128gb:ncpus=10
#PBS -l walltime=48:00:00
#PBS -N macs2_IgG_heatmap
#PBS -J 0-11

tid=$PBS_ARRAY_INDEX

module load anaconda3/personal
source activate deeptools

DIR=/rds/general/user/mml120/ephemeral/p2/CR_000
OUTDIR=/rds/general/user/mml120/ephemeral/p2/CR_000/heatmap/macs2/unnorm_IgG/narrowPeak

files=$(ls $DIR/alignment/bigwig/CR_000_*_raw.bw)
arr=($files)
sample=${arr[$tid]}

sampleID=$(echo `basename $sample _raw.bw`)

cores=8
computeMatrix reference-point -S $DIR/alignment/bigwig/"$sampleID"_raw.bw \
              -R $DIR/peakCalling/macs2/unnorm_Ig/"$sampleID"_peaks.narrowPeak \
              --skipZeros -o $OUTDIR/"$sampleID"_macs2_narrowPeak.mat.gz -p $cores -a 3000 -b 3000 --referencePoint center

plotHeatmap -m $OUTDIR/"$sampleID"_macs2_narrowPeak.mat.gz -out $OUTDIR/"$sampleID"_macs2_narrowPeak_heatmap.png --sortUsing sum --startLabel "Peak Start" -\
-endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "$sampleID"
