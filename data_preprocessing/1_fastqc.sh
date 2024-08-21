#!/bin/bash

#PBS -l select=1:mem=64gb:ncpus=1
#PBS -l walltime=04:00:00
#PBS -N fastqc_CR_000
#PBS -J 0-31

tid=$PBS_ARRAY_INDEX

DIR=/rds/general/user/mml120/ephemeral/p2/CR_000/

module load anaconda3/personal

source activate cutnTag

files=$(ls $DIR/raw/X204SC24022143-Z01-F001/01.RawData/*/*.fq.gz)

arr=($files)

sample=${arr[$tid]}

fastqc $sample --dir . -o .

mv * $DIR/fastqc
