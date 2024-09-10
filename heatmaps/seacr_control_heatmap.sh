##seacr_control_heatmap.sh is a script used to generate heatmaps for SEACR CUT&Run control peaks

#!/bin/bash
#PBS -l select=1:mem=128gb:ncpus=10
#PBS -l walltime=48:00:00
#PBS -N peak_heatmap_CR_000
#PBS -J 0-11

tid=$PBS_ARRAY_INDEX

module load anaconda3/personal
source activate deeptools

DIR=/path/to/home/directory/alignment/bigwig
peaks=/path/to/home/directory/peakCalling
OUTDIR=/path/to/home/directory/heatmap/seacr

files=$(ls $DIR/CR_000_*_raw.bw)
arr=($files)
sample=${arr[$tid]}

sampleID=$(echo `basename $sample _raw.bw`)

awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' $peaks/SEACR/"$sampleID"_seacr_control.peaks.stringent.bed >$peaks/SEACR_1/"$sampleID"_seacr_control.peaks.summitRegion.bed

cores=8
computeMatrix reference-point -S $DIR/"$sampleID"_raw.bw \
              -R $peaks/SEACR_1/"$sampleID"_seacr_control.peaks.summitRegion.bed \
              --skipZeros -o $OUTDIR/"$sampleID"_SEACR.mat.gz -p $cores -a 3000 -b 3000 --referencePoint center

plotHeatmap -m $OUTDIR/"$sampleID"_SEACR.mat.gz -out $OUTDIR/"$sampleID"_SEACR_heatmap.png --sortUsing sum --startLabel "Peak Start" -\
-endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "$sampleID"
