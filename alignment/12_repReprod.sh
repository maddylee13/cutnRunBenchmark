##12_repReprod.sh is a script used to study the reproducibility between replicates and across conditions

#!/bin/bash

#PBS -l select=1:mem=32gb:ncpus=8 
#PBS -l walltime=48:00:00
#PBS -N CR_000
#PBS -J 0-15

tid=$PBS_ARRAY_INDEX

module load anaconda3/personal
source activate cutnTag

projPath=/rds/general/user/mml120/ephemeral/p2/CR_000

files=$(ls $projPath/alignment/bed/*_bowtie2.fragments.bed)
arr=($files)
sample=${arr[$tid]}

sampleID=$(echo `basename $sample _bowtie2.fragments.bed`)

binLen=500

awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $projPath/alignment/bed/"$sampleID"_bowtie2.fragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  >$projPath/alignment/bed/"$sampleID"_bowtie2.fragmentsCount.bin$binLen.bed

