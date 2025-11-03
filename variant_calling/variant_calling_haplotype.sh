## VARIANT CALLING W/ HAPLOTYPECALLER

TUM_ID=$1; 
echo ${TUM_ID};

cd /stereoseq/all_samples/vcf/${TUM_ID}/;

/home/ubuntu/tools/gatk-4.4.0.0/gatk HaplotypeCaller \
-R /stereoseq/reference/GRCh38.109.fa \
-I /stereoseq/all_samples/bams/${TUM_ID}/${TUM_ID}.recal.RGA.bam \
-O ${TUM_ID}.all.vcf.gz \
--tmp-dir /stereoseq/tmp \
--dont-use-soft-clipped-bases \
2>${TUM_ID}.all.vcf.gz.log ; 
gzip -dkf ${TUM_ID}.all.vcf.gz; # unzip
# aws s3 cp ${TUM_ID}.all.vcf.gz s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/vcf_updatedFeb25/;
# aws s3 cp ${TUM_ID}.all.vcf.gz.log s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/vcf_updatedFeb25/;
# aws s3 cp ${TUM_ID}.all.vcf s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/vcf_updatedFeb25/;
