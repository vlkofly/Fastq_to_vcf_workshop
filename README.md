1. login to Metacentrum
ssh vlkofly@zuphux.cerit-sc.cz
2. create a workshop directory
mkdir fastq_to_vcf_workshop
cd fastq_to_vcf_workshop
3. get the data
cp -r /storage/brno3-cerit/home/vlkofly/mockingbird_genomics/fastq_vcf_workshop .

4. map the reads
python3 ../Fastq-to-vcf/1_mapreads.py -samplenames ../sample.list.mapping.txt -wd `pwd` -datadir ../fastq_data/ -fastqc t -trial t -print true -ref /storage/brno3-cerit/home/filip_kolar/JIC_reference/alygenomes.fasta

5. Quality check
what is depth of coverage?
percentage of mapped reads?

in logstat folder:
python36-modules-gcc	
multiqc .

5.  Run Haplotype caller

python3 ../../Fastq-to-vcf/2_callvars.py -ploidyfile ../../ploidy.file.txt -trial f -print false -ref /storage/brno3-cerit/home/filip_kolar/JIC_reference/alygenomes.fasta -workdir `pwd`

6. Run Joint genotyping

qsub -v 'samples=sample.list.txt,sourcedirs="/storage/brno3-cerit/home/vlkofly/mockingbird_genomics/fastq_vcf_workshop/workdir/HC",outdir=vcf,outvcf=arenosa.vcf.gz,NV=no' -N variant_call ../Fastq-to-vcf/3_genotypeGVCF.sh

7.  Run Filtering

qsub -v 'vcf_var=arenosa.vcf.gz' -N filter.vcf ../Fastq-to-vcf/4_filter.sh

