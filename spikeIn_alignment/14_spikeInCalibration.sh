##14_spikeInCalibration.sh is a script used to calibrate the data to the spike in genome

#!/bin/bash
 
#PBS -l select=1:mem=32gb:ncpus=8 
#PBS -l walltime=48:00:00
#PBS -N CR_000_spikeIn_calibration
#PBS -J 0-15
 
tid=$PBS_ARRAY_INDEX
 
module load anaconda3/personal
source activate cutnTag
 
projPath=/path/to/home/directory
GENOMESIZE=/path/to/genome/directory/chm13v2.0/chm13v2.0.chrom.sizes

files=$(ls $projPath/alignment/bed/*_bowtie2.fragments.bed)
arr=($files)
sample=${arr[$tid]}

files2=$(ls $projPath/alignment/sam/*_bowtie2_spikeIn.sam)
arr2=($files2)
sample2=${arr2[$tid2]}

sampleID=$(echo `basename $sample _bowtie2.fragments.bed`)

sampleID2=$(echo `basename $sample2 _bowtie2_spikeIn.sam`)

seqDepthDouble=`samtools view -F 0x04 $projPath/alignment/sam/"$sampleID2"_bowtie2_spikeIn.sam | wc -l`
seqDepth=$((seqDepthDouble/2))


if [[ "$seqDepth" -gt "1" ]]; then

mkdir -p $projPath/alignment/bedgraph

scale_factor=`echo "10000 / $seqDepth" | bc -l`
echo "Scaling factor for $sampleID is: $scale_factor!"
bedtools genomecov -bg -scale $scale_factor -i $projPath/alignment/bed/"$sampleID"_bowtie2.fragments.bed -g $GENOMESIZE> $projPath/alignment/bedgraph/"$sampleID"_bowtie2.fragments.normalized.bedgraph

fi

