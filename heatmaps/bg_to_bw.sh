##bg_to_bw.sh is a script to convert .bedgraph files to .bw

#!/bin/bash
#PBS -l select=1:mem=64gb:ncpus=1
#PBS -l walltime=48:00:00
#PBS -N bigwig_norm_CR_000
#PBS -J 0-15

tid=$PBS_ARRAY_INDEX

module load anaconda3/personal
source activate ucsc-tools

DIR=/path/to/home/directory/alignment/bedgraph/normalized
OUTDIR=/path/to/home/directory/peakCalling/normalized/bigwig
chromSizes=/path/to/genome/directory/genome_indexes/chm13v2.0/chm13v2.0.chrom.sizes

files=$(ls $DIR/*_bowtie2.fragments.normalized.bedgraph)
arr=($files)
sample=${arr[$tid]}

sampleID=$(echo `basename $sample _bowtie2.fragments.normalized.bedgraph`)

bedGraphToBigWig $DIR/"$sampleID"_bowtie2.fragments.normalized.bedgraph $chromSizes $OUTDIR/"$sampleID"_normalized.bw
