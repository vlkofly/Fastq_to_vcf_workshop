#!/bin/bash
#PBS -N  p2m1m2_HRA_01ta
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=3:mem=20gb:scratch_local=5gb
#PBS -j oe

trap 'clean_scratch' TERM EXIT

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 
cp /storage/brno3-cerit/home/vlkofly/mockingbird_genomics/fastq_vcf_workshop/mapping/../fastq_data/Arenosa2_1000L_R1.fastq.gz $SCRATCHDIR
cp /storage/brno3-cerit/home/vlkofly/mockingbird_genomics/fastq_vcf_workshop/mapping/../fastq_data/Arenosa2_1000L_R2.fastq.gz $SCRATCHDIR
cp -r /storage/brno3-cerit/home/filip_kolar/JIC_reference $SCRATCHDIR/
echo input files copied at `date`

cd $SCRATCHDIR
module add trimmomatic-0.36
module add bwa-0.7.15
module add samtools-1.4
module add picard-2.8.1 
module add parallel-20160622
module add fastQC-0.11.5 
module add parallel-20160622

mkdir outtrim
mkdir logstat

cp /storage/brno3-cerit/home/filip_kolar/600_genomes/adapters/TruSeq3-PE-2.fa $SCRATCHDIR

java -XX:ParallelGCThreads=3 -Xmx20g -jar /software/trimmomatic/0.36/dist/jar/trimmomatic-0.36.jar PE -threads 2 -trimlog logstat/trim_lane0Arenosa2.log Arenosa2_1000L_R1.fastq.gz Arenosa2_1000L_R2.fastq.gz -baseout trim_Arenosa2_1000L_R1.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:23:10 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:50

mv trim_* outtrim/
echo trimming of lane 0 finished at `date`

echo "processing lane 0"  
lanename=`zcat Arenosa2_1000L_R1.fastq.gz | head -n 1 | cut -d: -f1,2,3,4`
bwa mem -M -R '@RG\tID:'"$lanename"'\tLB:'"$lanename"'one\tPL:illumina\tSM:HRA_01ta' -t 3 JIC_reference/alygenomes.fasta outtrim/trim_*Arenosa2_1000L_R1*_1P.fastq.gz outtrim/trim_*Arenosa2_1000L_R1*_2P.fastq.gz | samtools view -@ 3 -bu - | samtools sort - -@ 3 -m 5G -o sort_Arenosa2_1000L_R1_paired.bam || exit 1  
echo paired reads mapped at `date`

bwa mem -M -R '@RG\tID:'"$lanename"'\tLB:'"$lanename"'one\tPL:illumina\tSM:HRA_01ta' -t 3 JIC_reference/alygenomes.fasta outtrim/trim_*Arenosa2_1000L_R1*_1U.fastq.gz | samtools view -@ 3 -bu - | samtools sort - -@ 3 -m 5G -o sort_Arenosa2_1000L_R1_unpaired1.bam || exit 1 
echo unpaired reads1 mapped at `date`

bwa mem -M -R '@RG\tID:'"$lanename"'\tLB:'"$lanename"'one\tPL:illumina\tSM:HRA_01ta' -t 3 JIC_reference/alygenomes.fasta outtrim/trim_*Arenosa2_1000L_R1*_2U.fastq.gz | samtools view -@ 3 -bu - | samtools sort - -@ 3 -m 5G -o sort_Arenosa2_1000L_R1_unpaired2.bam || exit 1 
echo unpaired reads2 lane $l mapped at `date`

samtools index sort_Arenosa2_1000L_R1_paired.bam
samtools index sort_Arenosa2_1000L_R1_unpaired1.bam
samtools index sort_Arenosa2_1000L_R1_unpaired2.bam
echo all reads of lane $l sorted `date`


for b in *bam
do
blist+="I=$b "
 done 
 echo $blist

java -XX:ParallelGCThreads=3 -Xmx20g -jar $PICARD281 MergeSamFiles $blist O=Arenosa2merged.bam 
java -XX:ParallelGCThreads=3 -Xmx20g -jar $PICARD281 BuildBamIndex I=Arenosa2merged.bam
echo all reads sorted and merged `date`

limit=$(echo `ulimit -n` - 50 | bc)
java -XX:ParallelGCThreads=3 -Xmx20g -jar $PICARD281 MarkDuplicates I=Arenosa2merged.bam O=Arenosa2dedup.bam MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=$limit M=logstat/Dup_metricsArenosa2.log ASSUME_SORTED=true TAGGING_POLICY=All
echo deduplicated `date`

java -XX:ParallelGCThreads=3 -Xmx20g -jar $PICARD281 BuildBamIndex I=Arenosa2dedup.bam
echo all reads sorted and merged `date`

parallel -j 2 "samtools {} Arenosa2dedup.bam > logstat/Arenosa2.{}.txt" ::: idxstats flagstat stats

samtools depth Arenosa2dedup.bam | awk '{sum+=$3} END {print sum/NR}' > logstat/depth_Arenosa2.txt 
samtools depth -b JIC_reference/LyV2_genes_only.bed  Arenosa2dedup.bam | awk '{sum+=$3} END {print sum/NR}' > logstat/codingdepth_Arenosa2.txt 
echo stats done at `date`

mv Arenosa2dedup.bam HRA_01ta_final.bam
mv Arenosa2dedup.bai HRA_01ta_final.bai

rm -r `dirname JIC_reference/alygenomes.fasta` 
rm *paired*.bam*
rm *.fastq.gz
rm *merged.ba*

cp -r * /storage/brno3-cerit/home/vlkofly/mockingbird_genomics/fastq_vcf_workshop/mapping/p2m1m2result/ || export CLEAN_SCRATCH=false
echo results exported at `date`
