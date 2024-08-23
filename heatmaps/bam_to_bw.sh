##bam_to_bw.sh is a script used to convert .bam files to .bw files

#!/bin/bash
#PBS -l select=1:mem=128gb:ncpus=10
#PBS -l walltime=48:00:00
#PBS -N bam_to_bigwig_CR_000
#PBS -J 0-15

tid=$PBS_ARRAY_INDEX

module load anaconda3/personal
source activate deeptools

DIR=/rds/general/user/mml120/ephemeral/p2/CR_000/
OUTDIR=/rds/general/user/mml120/ephemeral/p2/CR_000/alignment/bigwig

files=$(ls $DIR/alignment/bam/*_bowtie2.mapped.bam)
arr=($files)
sample=${arr[$tid]}

sampleID=$(echo `basename $sample _bowtie2.mapped.bam`)

samtools sort -o $DIR/alignment/bam/"$sampleID"_bowtie2.sorted.mapped.bam $DIR/alignment/bam/"$sampleID"_bowtie2.mapped.bam
samtools index $DIR/alignment/bam/"$sampleID"_bowtie2.sorted.mapped.bam

bamCoverage -b $DIR/alignment/bam/"$sampleID"_bowtie2.sorted.mapped.bam -o $OUTDIR/"$sampleID"_raw.bw
