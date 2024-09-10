##10_fileFormat.sh is a script used to convert files, that haven't had duplicates removed, from .sam to .bam, as well as .bam to .bed

#!/bin/bash
#PBS -l select=1:mem=100gb:ncpus=8 
#PBS -l walltime=48:00:00
#PBS -N CR_000_fileFormat_not_remdup
#PBS -J 0-15

tid=$PBS_ARRAY_INDEX

module load anaconda3/personal
source activate cutnTag

DIR=/path/to/home/directory/alignment

files=$(ls $DIR/sam/bowtie2_summary/*_bowtie2.sam)
arr=($files)
sample=${arr[$tid]}

sampleID=$(echo `basename $sample _bowtie2.sam`)

samtools view -bS -F 0x04 "$DIR/sam/bowtie2_summary/"$sampleID"_bowtie2.sam" > "$DIR/bam/"$sampleID"_bowtie2.mapped.bam"

bedtools bamtobed -i $DIR/bam/"$sampleID"_bowtie2.mapped.bam -bedpe > "$DIR/bed/"$sampleID"_bowtie2.bed"

awk '$1==$4 && $6-$2 < 1000 {print $0}' $DIR/bed/"$sampleID"_bowtie2.bed >$DIR/bed/"$sampleID"_bowtie2.clean.bed

cut -f 1,2,6 $DIR/bed/"$sampleID"_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$DIR/bed/"$sampleID"_bowtie2.fragments.bed
