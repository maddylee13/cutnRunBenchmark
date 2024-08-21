##4_spikeIn_bowtie2.sh is a script used to align the data to the E. coli genome

#!/bin/bash
 
#PBS -l select=1:mem=64gb:ncpus=10
#PBS -l walltime=48:00:00
#PBS -N CR_000_eAlign
#PBS -J 0-15

tid=$PBS_ARRAY_INDEX

module load anaconda3/personal
source activate bowtie2

DIR=/rds/general/user/mml120/home/p2/fastp

GENOMEDIR=/rds/general/project/traditiom/live/epimem/cutntag/spikein/ecoli/Escherichia_coli_K_12_MG1655/NCBI/2001-10-15/Sequence/Bowtie2Index

TEMPDIR=/rds/general/user/mml120/home/p2 

files=$(ls $DIR/*_1.trimmed.fastq.gz)
arr=($files)
sample=${arr[$tid]}

sampleID=$(echo `basename $sample _1.trimmed.fastq.gz`)

in1=$DIR/"$sampleID"_1.trimmed.fastq.gz
in2=$DIR/"$sampleID"_2.trimmed.fastq.gz

cores=8

bowtie2 --local --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 -I 10 -X 700 -p $cores -x $GENOMEDIR/genome -1 $in1 -2 $in2 -S $TEMPDIR/alignment/sam/"$sampleID"_bowtie2_spikeIn.sam &> $TEMPDIR/alignment/sam/bowtie2_summary/"$sampleID"_bowtie2_spikeIn.txt

seqDepthDouble=$(samtools view -F 0x04 -c $TEMPDIR/alignment/sam/"$sampleID"_bowtie2_spikeIn.sam)
seqDepth=$((seqDepthDouble / 2))
echo $seqDepth >$TEMPDIR/alignment/sam/bowtie2_summary/"$sampleID"_bowtie2_spikeIn.seqDepth

