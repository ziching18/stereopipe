# BULK RNASEQ PRE-PROCESSING (FQ)

TUM_ID=$1; 

### trim
cd /stereoseq/all_samples/rnaseq/$TUM_ID/;
trim_galore ${TUM_ID}_r1.fq.gz ${TUM_ID}_r2.fq.gz --paired --gzip -o . -a CTGTCTCTTATACACATCT 2>$TUM_ID.log ; 
aws s3 cp ${TUM_ID}_r1_val_1.fq.gz s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/1_trim_galore_out/
aws s3 cp ${TUM_ID}_r2_val_2.fq.gz s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/1_trim_galore_out/
rm ${TUM_ID}._r1_val_1.fq.gz ; 
rm ${TUM_ID}._r2_val_2.fq.gz ; 


### map
cd /stereoseq/all_samples/rnaseq/$TUM_ID/;
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/1_trim_galore_out/${TUM_ID}_r1_val_1.fq.gz ./;
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/1_trim_galore_out/${TUM_ID}_r2_val_2.fq.gz ./;
zcat ${TUM_ID}_r1_val_1.fq.gz | awk '(NR%2==0){$0=substr($0,1,75)}{print}' > ${TUM_ID}.r_1_bwa_in.fq ; 
zcat ${TUM_ID}_r2_val_2.fq.gz | awk '(NR%2==0){$0=substr($0,1,75)}{print}' > ${TUM_ID}.r_2_bwa_in.fq ; 
bwa mem -M -t 4 /stereoseq/reference/GRCh38.109.bwa.fa \
${TUM_ID}.r_1_bwa_in.fq \
${TUM_ID}.r_2_bwa_in.fq 2>${TUM_ID}.bwa.log | \
samtools view --threads 8 -b - | \
samtools sort --threads 8 > ${TUM_ID}.RNAseq.sorted.bam ; 
aws s3 cp ${TUM_ID}.bwa.log s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/2_bwa_out/ ; \
aws s3 cp ${TUM_ID}.RNAseq.sorted.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/2_bwa_out/ && touch ${TUM_ID}.map.success || touch ${TUM_ID}.map.failed ; 
aws s3 cp ${TUM_ID}.r_1_bwa_in.fq s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/2_bwa_out/ ; \
aws s3 cp ${TUM_ID}.r_2_bwa_in.fq s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/2_bwa_out/; 
rm *sorted.bam;
rm ${TUM_ID}.r_1_bwa_in.fq  ;
rm ${TUM_ID}.r_2_bwa_in.fq  ;
rm ${TUM_ID}._r1_val_1.fq.gz ; 
rm ${TUM_ID}._r2_val_2.fq.gz ; 
rm cp ${TUM_ID}.RNAseq.sorted.bam ;
cp ${TUM_ID}.bwa.log ;


### add read groups
cd /stereoseq/all_samples/rnaseq/$TUM_ID/;
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/2_bwa_out/${TUM_ID}.RNAseq.sorted.bam ./;
/home/ubuntu/tools/gatk-4.4.0.0/gatk AddOrReplaceReadGroups \
-I ${TUM_ID}.RNAseq.sorted.bam \
-O ${TUM_ID}.RNAseq.sorted.RGA.bam \
--TMP_DIR /stereoseq/tmp \
--RGID ${TUM_ID} \
--RGLB RNAseq \
--RGPL Illumina \
--RGPU RNAseq.${TUM_ID} \
--RGSM ${TUM_ID}_RNAseq_TUM \
&& touch ${TUM_ID}.addr.success || touch ${TUM_ID}.addr.failed ; 
aws s3 cp ${TUM_ID}.RNAseq.sorted.RGA.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/3_rga_out/ ; 
rm ${TUM_ID}.RNAseq.sorted.RGA.bam;
rm ${TUM_ID}.RNAseq.sorted.bam;


## mark duplicates
cd /stereoseq/all_samples/rnaseq/$TUM_ID/;
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/3_rga_out/${TUM_ID}.RNAseq.sorted.RGA.bam ./;
/home/ubuntu/tools/gatk-4.4.0.0/gatk MarkDuplicates \
-I ${TUM_ID}.RNAseq.sorted.RGA.bam \
-O ${TUM_ID}.RNAseq.dedup.bam \
-M ${TUM_ID}.MarkDup.metrics \
--CREATE_INDEX true \
--VALIDATION_STRINGENCY SILENT ; 
aws s3 cp ${TUM_ID}.RNAseq.dedup.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/4_dedupped_bams/ ; 
aws s3 cp ${TUM_ID}.RNAseq.dedup.bai s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/4_dedupped_bams/ ; 
aws s3 cp ${TUM_ID}.MarkDup.metrics s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/4_dedupped_bams/ ; 
rm ${TUM_ID}.RNAseq.dedup.bam;
rm ${TUM_ID}.RNAseq.dedup.bai;
rm ${TUM_ID}.MarkDup.metrics;
rm ${TUM_ID}.RNAseq.sorted.RGA.bam;


## split reads
cd /stereoseq/all_samples/rnaseq/$TUM_ID/;
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/4_dedupped_bams/${TUM_ID}.RNAseq.dedup.bam ./ ; 
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/4_dedupped_bams/${TUM_ID}.RNAseq.dedup.bai ./ ; 
/home/ubuntu/tools/gatk-4.4.0.0/gatk SplitNCigarReads \
-I ${TUM_ID}.RNAseq.dedup.bam \
-O ${TUM_ID}.RNAseq.split.bam \
-R /stereoseq/reference/GRCh38.109.fa \
--tmp-dir /stereoseq/tmp \
&& touch ${TUM_ID}.split.success || touch ${TUM_ID}.split.failed ; 
aws s3 cp ${TUM_ID}.RNAseq.split.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/5_split_bams/ ; 
aws s3 cp ${TUM_ID}.RNAseq.split.bai s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/5_split_bams/ ; 
# rm ${TUM_ID}.RNAseq.split.bam;
# rm ${TUM_ID}.RNAseq.split.bai;
rm ${TUM_ID}.RNAseq.dedup.bam;
rm ${TUM_ID}.RNAseq.dedup.bai;


## base quality recalibration - split
cd /stereoseq/all_samples/rnaseq/$TUM_ID/;
# aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/5_split_bams/${TUM_ID}.RNAseq.split.bam ./ ; 
# aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/5_split_bams/${TUM_ID}.RNAseq.split.bai ./ ; 
/home/ubuntu/tools/gatk-4.4.0.0/gatk BaseRecalibrator \
-I ${TUM_ID}.RNAseq.split.bam \
-O ${TUM_ID}.RNAseq.recal_data.grp \
-R /stereoseq/reference/GRCh38.109.fa \
--tmp-dir /stereoseq/tmp \
--known-sites /stereoseq/reference/Homo_sapiens_assembly38.known_indels.renamed.vcf \
--known-sites /stereoseq/reference/Homo_sapiens_assembly38.dbsnp138.renamed.vcf \
--known-sites /stereoseq/reference/Mills_and_1000G_gold_standard.indels.hg38.renamed.vcf \
&& touch ${TUM_ID}.recal_1.success || touch ${TUM_ID}.recal_1.failed ; 
aws s3 cp ${TUM_ID}.RNAseq.recal_data.grp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/5_split_bams/ ; 
# rm ${TUM_ID}.RNAseq.split.bam;
# rm ${TUM_ID}.RNAseq.split.bai;


## base quality recalibration - apply BQSR
cd /stereoseq/all_samples/rnaseq/$TUM_ID/;
# aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/5_split_bams/${TUM_ID}.RNAseq.split.bam ./ ; 
# aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/5_split_bams/${TUM_ID}.RNAseq.split.bai ./ ; 
/home/ubuntu/tools/gatk-4.4.0.0/gatk ApplyBQSR \
-I ${TUM_ID}.RNAseq.split.bam \
-O ${TUM_ID}.RNAseq.recal.bam \
-R /stereoseq/reference/GRCh38.109.fa \
-bqsr ${TUM_ID}.RNAseq.recal_data.grp \
--tmp-dir /stereoseq/tmp \
&& touch ${TUM_ID}.recal_2.success || touch ${TUM_ID}.recal_2.failed ; 
aws s3 cp ${TUM_ID}.RNAseq.recal.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/5_split_bams/ ; 
aws s3 cp ${TUM_ID}.RNAseq.recal.bai s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/rnaseq_tum/5_split_bams/ ; 
# rm ${TUM_ID}.RNAseq.recal.bam;
# rm ${TUM_ID}.RNAseq.recal.bai;

