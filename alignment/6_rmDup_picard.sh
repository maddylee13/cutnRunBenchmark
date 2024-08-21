##6_rmDup_picard.sh is a script used to remove duplicated fragments

#!/bin/bash

#PBS -l select=1:mem=64gb:ncpus=10
#PBS -l walltime=48:00:00
#PBS -N _remdup_picard_CR_000
#PBS -J 0-15

tid=$PBS_ARRAY_INDEX

module load anaconda3/personal
source activate cutnTag


DIR=/rds/general/user/mml120/ephemeral/p2/CR_000/alignment

PICARD=/rds/general/user/mml120/home/anaconda3/pkgs/picard-2.18.7-2/share/picard-2.18.7-2/picard.jar

ephem=/rds/general/user/mml120/ephemeral

files=$(ls $DIR/sam/bowtie2_summary/*_bowtie2.sam)
arr=($files)
sample=${arr[$tid]}

sampleID=$(echo `basename $sample _bowtie2.sam`)

java -jar $PICARD SortSam I=$DIR/sam/bowtie2_summary/"$sampleID"_bowtie2.sam O=$DIR/sam/"$sampleID"_bowtie2.sorted.sam SORT_ORDER=coordinate

java -jar $PICARD MarkDuplicates I=$DIR/sam/"$sampleID"_bowtie2.sorted.sam O=$DIR/removeDuplicate/"$sampleID"_bowtie2.sorted.dupMarked.sam METRICS_FILE=$DIR/removeDuplicate/picard_summary/"$sampleID"_picard.dupMark.txt TMP_DIR=$ephem/TMP

java -jar $PICARD MarkDuplicates I=$DIR/sam/"$sampleID"_bowtie2.sorted.sam O=$DIR/removeDuplicate/"$sampleID"_bowtie2.sorted.rmDup.sam REMOVE_DUPLICATES=true METRICS_FILE=$DIR/removeDuplicate/picard_summary/"$sampleID"_picard.rmDup.txt TMP_DIR=$ephem/TMP

#checking files

cd $DIR/removeDuplicate

module load anaconda3/personal

samtools quickcheck -v *_bowtie2.sorted.rmDup.sam > bad_sams.fofn   && echo 'all ok' || echo 'some files failed check, see bad_sams.fofn'

ls -lh $DIR/sam/*.sorted.sam 
ls -lh $DIR/removeDuplicate/*.sorted.dupMarked.sam 
ls -lh $DIR/removeDuplicate/picard_summary/*_picard.dupMark.txt 
ls -lh $DIR/removeDuplicate/*_bowtie2.sorted.rmDup.sam 
