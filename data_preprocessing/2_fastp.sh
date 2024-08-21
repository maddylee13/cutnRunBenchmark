##2_fastp.sh is a script used to trim the adaptors using the fastp programme

#!/bin/bash

#PBS -l select=1:mem=30gb:ncpus=6
#PBS -l walltime=02:00:00
#PBS -N fastp_CR_000
#PBS -J 0-31
tid=$PBS_ARRAY_INDEX

DIR=/rds/general/user/mml120/ephemeral/p2/CR_000

OUTDIR=/rds/general/user/mml120/ephemeral/p2/CR_000/fastp

module load anaconda3/personal
source activate cutnTag

files=$(ls $DIR/*/*/*/*/*.fq.gz)
arr=($files)
sample=${arr[$tid]} 

sampleID=$(echo `basename $sample _1.fq.gz`)

WORKDIR=$(echo `dirname $sample`)

fastp -i $DIR/*/*/*/*/*"$sampleID"_1.fq.gz -I $DIR/*/*/*/*/*"$sampleID"_2.fq.gz -o "$sampleID"_1.trimmed.fastq.gz -O "$sampleID"_2.trimmed.fastq.gz --detect_adapter_for_pe -l 20 -j "$sampleID".fastp.json -h "$sampleID".fastp.html

mv * $OUTDIR
