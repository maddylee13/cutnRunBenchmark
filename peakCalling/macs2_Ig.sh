##macs2_Ig.sh is a script used to call peaks using the MACS2 peak caller with an IgG control sample

#!/bin/bash

#PBS -l select=1:mem=24gb:ncpus=2
#PBS -l walltime=48:00:00
#PBS -N CR_000_macs2_Ig
#PBS -J 0-11

module load anaconda3/personal
source activate macs2

tid=$PBS_ARRAY_INDEX

cd $DIR/peakCalling/macs2

DIR=/path/to/home/directory

for file in $DIR/alignment/bam/noIgG/*_bowtie2.mapped.bam
do
  sampleID=$(basename $file _bowtie2.mapped.bam)
  OUTDIR=$DIR/peakCalling/macs2/unnorm_Ig # Define the output directory
  log=$OUTDIR/logs
  IgG=/rds/general/user/mml120/ephemeral/p2/CR_000/alignment/bam/IgG/CR_000_I2A_bowtie2.mapped.bam

  macs2 callpeak -t $file -c $IgG -f BAMPE --keep-dup all --outdir $OUTDIR 2> "${log}/${sampleID}_macs2.log" -n $sampleID  --nomodel --shift -75 --extsize 150 
done

