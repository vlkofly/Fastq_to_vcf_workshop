#!/bin/bash 
#PBS -N V1.HRA_10ta
#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=2:mem=16gb:scratch_local=4gb
#PBS -j oe

module add gatk-3.7 
module add samtools-1.4 
trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 
cp -r /storage/brno3-cerit/home/filip_kolar/JIC_reference $SCRATCHDIR
cp /storage/brno3-cerit/home/vlkofly/mockingbird_genomics/fastq_vcf_workshop/workdir/p2m1m2result/HRA_10ta*final.ba* $SCRATCHDIR
cd $SCRATCHDIR
echo data loaded at `date`

java -XX:ParallelGCThreads=2 -Xmx16g -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -I HRA_10ta_final.bam --min_base_quality_score 20 --min_mapping_quality_score 20 -rf BadMate -R JIC_reference/alygenomes.fasta -o HRA_10ta_HC.g.vcf.gz -ploidy 4 -stand_call_conf 20 -ERC GVCF --pcr_indel_model NONE -nct 2 --max_num_PL_values 350 
echo calling done at `date`

rm -r JIC_reference/alygenomes.fasta
rm HRA_10ta*final.ba* 
cp *vcf* /storage/brno3-cerit/home/vlkofly/mockingbird_genomics/fastq_vcf_workshop/workdir/p2m1m2result/../HC/ || export CLEAN_SCRATCH=false
printf "\nFinished\n\n"
