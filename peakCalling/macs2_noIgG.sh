##macs2_noIgG.sh is a script that calls peaks using the MACS2 peak caller without IgG as a control

#PBS -l select=1:mem=24gb:ncpus=2
#PBS -l walltime=48:00:00
#PBS -N CR_000_macs2
#PBS -J 0-15

module load anaconda3/personal
source activate macs2

tid=$PBS_ARRAY_INDEX

DIR=/path/to/home/directory

cd $DIR/peakCalling/macs2

for file in $DIR/alignment/bam/*_bowtie2.mapped.bam
do
  sampleID=$(basename $file _bowtie2.mapped.bam)
  OUTDIR=$DIR/peakCalling/macs2/unnorm_noIg 
  log=$OUTDIR/logs

  macs2 callpeak -t $file -f BAMPE --keep-dup all --outdir $OUTDIR 2> "${log}/${sampleID}_macs2.log" -n $sampleID  --nomodel --shift -75 --extsize 150 
done
