##11_fileFormat_rmDup.sh is a script used to convert files, with removed duplicates, from .sam to .bam, as well as .bam to .bed

#PBS -l select=5:mem=124gb:ncpus=10
#PBS -l walltime=48:00:00
#PBS -N fileFormat_rmDup
#PBS -J 0-3

tid=$PBS_ARRAY_INDEX

module load anaconda3/personal
source activate cutnTag

DIR=/path/to/home/directory/alignment

files=$(ls $DIR/removeDuplicate/forward/*_bowtie2.sorted.rmDup.sam)
arr=($files)
sample=${arr[$tid]}

sampleID=$(echo `basename $sample _bowtie2.sorted.rmDup.sam`)

samtools view -bS -F 0x04 $DIR/removeDuplicate/forward/"$sampleID"_bowtie2.sorted.rmDup.sam >$DIR/bam/"$sampleID"_bowtie2.sorted.rmDup.mapped.bam

bedtools bamtobed -i $DIR/bam/"$sampleID"_bowtie2.sorted.rmDup.mapped.bam -bedpe >$DIR/bed/"$sampleID"_rmDup_bowtie2.bed

awk '$1==$4 && $6-$2 < 1000 {print $0}' $DIR/bed/"$sampleID"_rmDup_bowtie2.bed >$DIR/bed/"$sampleID"_rmDup_bowtie2.clean.bed

cut -f 1,2,6 $DIR/bed/"$sampleID"_rmDup_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$DIR/bed/"$sampleID"_rmDup_bowtie2.fragments.bed
