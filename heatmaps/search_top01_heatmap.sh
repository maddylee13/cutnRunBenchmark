##search_top01_heatmap.sh is a script used to generate heatmaps for SEACR CUT&Run top 1% peaks

#!/bin/bash
#PBS -l select=1:mem=128gb:ncpus=10
#PBS -l walltime=48:00:00
#PBS -N peak_heatmap_CR_000
#PBS -J 0-11

tid=$PBS_ARRAY_INDEX

module load anaconda3/personal
source activate deeptools

DIR=/rds/general/user/mml120/ephemeral/p2/CR_000/alignment/bigwig
peaks=/rds/general/user/mml120/ephemeral/p2/CR_000/peakCalling
OUTDIR=/rds/general/user/mml120/ephemeral/p2/CR_000/heatmap/seacr/top001

files=$(ls $DIR/CR_000_*_raw.bw)
arr=($files)
sample=${arr[$tid]}

sampleID=$(echo `basename $sample _raw.bw`)

awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' $peaks/SEACR_peaks/"$sampleID"_seacr_top0.01.peaks.stringent.bed >$peaks/SEACR_peaks/"$sampleID"_seacr_top0.01.peaks.summitRegion.bed

cores=8
computeMatrix reference-point -S $DIR/"$sampleID"_raw.bw \
              -R $peaks/SEACR_peaks/"$sampleID"_seacr_top0.01.peaks.summitRegion.bed \
              --skipZeros -o $OUTDIR/"$sampleID"_SEACR_top001.mat.gz -p $cores -a 3000 -b 3000 --referencePoint center

plotHeatmap -m $OUTDIR/"$sampleID"_SEACR_top001.mat.gz -out $OUTDIR/"$sampleID"_SEACR_top001_heatmap.png --sortUsing sum --startLabel "Peak Start" -\
-endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "$sampleID"
