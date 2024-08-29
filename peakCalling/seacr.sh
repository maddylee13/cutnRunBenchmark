##seacr.sh is a script used to call peaks using the SEACR peak caller

#!/bin/bash
 
#PBS -l select=1:mem=24gb:ncpus=2
#PBS -l walltime=48:00:00
#PBS -N CR_000_SEACR
#PBS -J 0-11

tid=$PBS_ARRAY_INDEX

DIR=/rds/general/user/mml120/ephemeral/p2/CR_000

module load anaconda3/personal
source activate seacr

files=$(ls $DIR/alignment/bedgraph/*_bowtie2.fragments.normalized.bedgraph)
arr=($files)
sample=${arr[$tid]}

sampleID=$(echo `basename $sample _bowtie2.fragments.normalized.bedgraph`)

seacr="/rds/general/user/mml120/home/anaconda3/pkgs/seacr-1.3-hdfd78af_2/bin/SEACR_1.3.sh"
histControl=$DIR/alignment/bedgraph/histControl/CR_000_I2A_bowtie2.fragments.normalized.bedgraph
mkdir -p $DIR/peakCalling/SEACR

bash $seacr $DIR/alignment/bedgraph/"$sampleID"_bowtie2.fragments.normalized.bedgraph \
     $histControl \
     non stringent $DIR/peakCalling/SEACR/"$sampleID"_seacr_control.peaks

bash $seacr $DIR/alignment/bedgraph/"$sampleID"_bowtie2.fragments.normalized.bedgraph 0.01 non stringent $DIR/peakCalling/SEACR/"$sampleID"_seacr_top0.01.peaks
