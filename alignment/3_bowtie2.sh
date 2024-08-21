##3_bowtie2.sh is a script used to align the data to the T2T genome

#!/bin/bash
 
#PBS -l select=1:mem=64gb:ncpus=10
#PBS -l walltime=48:00:00
#PBS -N  fastqc_CR_000_alignment
#PBS -J 0-15

DIR=/rds/general/user/mml120/home/p2/fastp
GENOMEDIR=/rds/general/project/traditiom/live/epimem/genome_indexes/chm13v2.0
TEMPDIR=/rds/general/user/mml120/home/p2

tid=$PBS_ARRAY_INDEX

module load anaconda3/personal
source activate bowtie2

files=$(ls $DIR/*_1.trimmed.fastq.gz)
arr=($files)
sample=${arr[$tid]}

sampleID=$(echo `basename $sample _1.trimmed.fastq.gz`)

in1=$DIR/"$sampleID"_1.trimmed.fastq.gz
in2=$DIR/"$sampleID"_2.trimmed.fastq.gz

cores=8

bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p $cores -x $GENOMEDIR/chm13v2.0 -1 $in1 -2 $in2 -S $TEMPDIR/alignment/sam/bowtie2_summary/"$sampleID"_bowtie2.sam &> $TEMPDIR/alignment/sam/bowtie2_summary/"$sampleID"_bowtie2.txt

