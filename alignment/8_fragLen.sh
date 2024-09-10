##8_fragLen.sh is a script to assess fragment size distribution

#!/bin/bash

#PBS -l select=5:mem=124gb:ncpus=10
#PBS -l walltime=48:00:00
#PBS -N fraglentgh_samtools
#PBS -J 0-15

tid=$PBS_ARRAY_INDEX

module load anaconda3/personal
source activate cutnTag

DIR=/path/to/home/directory/alignment

files=$(ls $DIR/sam/bowtie2_summary/*_bowtie2.sam)
arr=($files)
sample=${arr[$tid]}

sampleID=$(echo `basename $sample _bowtie2.sam`)

samtools view -F 0x04 $DIR/sam/bowtie2_summary/"$sampleID"_bowtie2.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >$DIR/sam/fragmentLen/"$sampleID"_fragmentLen.txt

