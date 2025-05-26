TUM_ID=$1;

cd /stereoseq/all_samples/vcf/${TUM_ID}/;

## VEP
/home/ubuntu/tools/ensembl-vep/vep \
--input_file /stereoseq/all_samples/vcf/${TUM_ID}/${TUM_ID}.somatic.funcotated.vcf \
--output_file /stereoseq/all_samples/vcf/${TUM_ID}/annotated_vcf/${TUM_ID}.somatic.funcotated.vep.vcf \
--format vcf --vcf --symbol --terms SO --tsl --biotype \
--hgvs --fasta /stereoseq/reference/GRCh38.109.fa \
--offline --cache --dir_cache /stereoseq/reference/ensembl_cache/ \
--plugin Frameshift --plugin Wildtype \
--dir_plugins /home/ubuntu/tools/VEP_plugins ;

## bgzip and index the pVACseq main input VCF
bgzip -c /stereoseq/all_samples/vcf/${TUM_ID}/annotated_vcf/${TUM_ID}.somatic.funcotated.vep.vcf \
> /stereoseq/all_samples/vcf/${TUM_ID}/annotated_vcf/${TUM_ID}.somatic.funcotated.vep.vcf.gz
tabix -p vcf /stereoseq/all_samples/vcf/${TUM_ID}/annotated_vcf/${TUM_ID}.somatic.funcotated.vep.vcf.gz ;

## Phased VCF
### Index all mutations VCF
tabix -p vcf /stereoseq/all_samples/vcf/${TUM_ID}/${TUM_ID}.all.vcf ;

### Sort all mutations VCF
/home/ubuntu/tools/gatk-4.4.0.0/gatk SortVcf \
-I /stereoseq/all_samples/vcf/${TUM_ID}/${TUM_ID}.all.vcf \
-O /stereoseq/all_samples/vcf/${TUM_ID}/${TUM_ID}.all.sorted.vcf \
-SD /stereoseq/reference/GRCh38.109.dict ;

### Create phased VCF w/ WhatsHap v2.3
whatshap phase \
/stereoseq/all_samples/vcf/${TUM_ID}/${TUM_ID}.all.sorted.vcf \
/stereoseq/all_samples/bams/${TUM_ID}/${TUM_ID}.recal.RGA.bam \
-o /stereoseq/all_samples/vcf/${TUM_ID}/annotated_vcf/${TUM_ID}.phased.vcf \
-r /stereoseq/reference/GRCh38.109.fa \
--output-read-list /stereoseq/all_samples/vcf/${TUM_ID}/annotated_vcf/${TUM_ID}.phased.read.list ;

### VEP annotate VCF
/home/ubuntu/tools/ensembl-vep/vep \
--input_file /stereoseq/all_samples/vcf/${TUM_ID}/annotated_vcf/${TUM_ID}.phased.vcf \
--output_file /stereoseq/all_samples/vcf/${TUM_ID}/annotated_vcf/${TUM_ID}.phased.vep.vcf \
--format vcf --vcf --symbol --terms SO --tsl \
--hgvs --fasta /stereoseq/reference/GRCh38.109.fa \
--offline --cache --dir_cache /stereoseq/reference/ensembl_cache/ \
--plugin Downstream --plugin Wildtype \
--dir_plugins /home/ubuntu/tools/VEP_plugins ;

### bgzip and index the phased VCF
# input into --phased-proximal-variants-vcf option
bgzip -c /stereoseq/all_samples/vcf/${TUM_ID}/annotated_vcf/${TUM_ID}.phased.vep.vcf \
> /stereoseq/all_samples/vcf/${TUM_ID}/annotated_vcf/${TUM_ID}.phased.vep.vcf.gz
tabix -p vcf /stereoseq/all_samples/vcf/${TUM_ID}/annotated_vcf/${TUM_ID}.phased.vep.vcf.gz 

## upload to AWS
# cd /stereoseq/all_samples/vcf/${TUM_ID}/annotated_vcf;
# aws s3 cp ${TUM_ID}.somatic.funcotated.vep.vcf s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/vcf_updatedFeb25/ ;
# aws s3 cp ${TUM_ID}.somatic.funcotated.vep.vcf.gz s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/vcf_updatedFeb25/ ;
# aws s3 cp ${TUM_ID}.phased.vcf s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/vcf_updatedFeb25/ ;
# aws s3 cp ${TUM_ID}.phased.vep.vcf s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/vcf_updatedFeb25/ ;
# aws s3 cp ${TUM_ID}.phased.vep.vcf.gz s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/vcf_updatedFeb25/ ;
# aws s3 cp ${TUM_ID}.phased.vep.vcf.gz.tbi s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/vcf_updatedFeb25/ ;
# aws s3 cp ${TUM_ID}.somatic.funcotated.vep.vcf.gz.tbi s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/vcf_updatedFeb25/ ;
