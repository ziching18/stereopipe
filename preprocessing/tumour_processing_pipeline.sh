## ALIGNED BAM FROM STEREOSEQ PRE-PROCESSING (BAM)

LIB_ID=$1; 
TUM_ID=$2; 
echo ${TUM_ID}

## merge bam files
aws s3 ls crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/1_sep_bams/ | 
grep ${TUM_ID}.bam$ | awk '{print $NF}'| while read l; do
    echo "-I "$l
done > ${TUM_ID}_bams.list

/home/ubuntu/tools/gatk-4.4.0.0/gatk MergeSamFiles --USE_THREADING true \
--arguments_file ${TUM_ID}_bams.list \
-O ${TUM_ID}.merged.bam

aws s3 cp ${TUM_ID}.merged.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/2_merged_bams/

cd /stereoseq/all_samples/bams/${TUM_ID}/;

## add read groups
echo `date "+%Y-%m-%d %H:%M:%S"` "AddOrReplaceReadGroups start" >> log.txt
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/2_merged_bams/${TUM_ID}.merged.bam ./ ; 
/home/ubuntu/tools/gatk-4.4.0.0/gatk AddOrReplaceReadGroups \
-I ${TUM_ID}.merged.bam \
-O ${TUM_ID}.merged.RGA.bam \
--TMP_DIR /stereoseq/tmp \
--RGID ${TUM_ID} \
--RGLB ${LIB_ID} \
--RGPL Illumina \
--RGPU ${LIB_ID}.${TUM_ID} \
--RGSM ${TUM_ID}_stereo_TUM \
&& touch ${TUM_ID}.addr.success || touch ${TUM_ID}.addr.failed ; 
aws s3 cp ${TUM_ID}.merged.RGA.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/2_merged_bams/ ; 
rm ${TUM_ID}.merged.bam; 
rm ${TUM_ID}.merged.RGA.bam;
echo `date "+%Y-%m-%d %H:%M:%S"` "AddOrReplaceReadGroups done" >> log.txt

## mark duplicates
echo `date "+%Y-%m-%d %H:%M:%S"` "MarkDuplicates start" >> log.txt
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/2_merged_bams/${TUM_ID}.merged.RGA.bam ./ ; 
/home/ubuntu/tools/gatk-4.4.0.0/gatk MarkDuplicates \
-I ${TUM_ID}.merged.RGA.bam \
-O ${TUM_ID}.dedup.bam \
-M ${TUM_ID}.MarkDup.metrics \
--TMP_DIR /stereoseq/tmp \
--CREATE_INDEX true \
--VALIDATION_STRINGENCY SILENT ; 
aws s3 cp ${TUM_ID}.dedup.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/3_dedupped_bams/ ; 
aws s3 cp ${TUM_ID}.dedup.bai s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/3_dedupped_bams/ ; 
aws s3 cp ${TUM_ID}.MarkDup.metrics s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/3_dedupped_bams/ ; 
rm ${TUM_ID}.merged.RGA.bam; 
rm ${TUM_ID}.dedup.bam; 
rm ${TUM_ID}.dedup.bai;
rm ${TUM_ID}.MarkDup.metrics;
echo `date "+%Y-%m-%d %H:%M:%S"` "MarkDuplicates done" >> log.txt


## split reads
echo `date "+%Y-%m-%d %H:%M:%S"` "SplitNCigarReads start" >> log.txt
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/3_dedupped_bams/${TUM_ID}.dedup.bam ./ ; 
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/3_dedupped_bams/${TUM_ID}.dedup.bai ./ ; 
/home/ubuntu/tools/gatk-4.4.0.0/gatk SplitNCigarReads \
-I ${TUM_ID}.dedup.bam \
-O ${TUM_ID}.split.bam \
-R /stereoseq/reference/GRCh38.109.fa \
--tmp-dir /stereoseq/tmp \
&& touch ${TUM_ID}.split.success || touch ${TUM_ID}.split.failed ; \
aws s3 cp ${TUM_ID}.split.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/4_split_bams/ ; 
aws s3 cp ${TUM_ID}.split.bai s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/4_split_bams/ ; 
rm ${TUM_ID}.dedup.bam; 
rm ${TUM_ID}.dedup.bai; 
rm ${TUM_ID}.split.bam; 
rm ${TUM_ID}.split.bai;
echo `date "+%Y-%m-%d %H:%M:%S"` "SplitNCigarReads done" >> log.txt

## base quality recalibration - split
echo `date "+%Y-%m-%d %H:%M:%S"` "BaseRecalibrator start" >> log.txt
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/4_split_bams/${TUM_ID}%.split.bam ./ ; 
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/4_split_bams/${TUM_ID}.split.bai ./ ; 
/home/ubuntu/tools/gatk-4.4.0.0/gatk BaseRecalibrator \
-I ${TUM_ID}.split.bam \
-O ${TUM_ID}.recal_data.grp \
-R /stereoseq/reference/GRCh38.109.fa \
--tmp-dir /stereoseq/tmp \
--known-sites /stereoseq/reference/Homo_sapiens_assembly38.known_indels.renamed.vcf \
--known-sites /stereoseq/reference/Homo_sapiens_assembly38.dbsnp138.renamed.vcf \
--known-sites /stereoseq/reference/Mills_and_1000G_gold_standard.indels.hg38.renamed.vcf \
&& touch ${TUM_ID}.recal_1.success || touch ${TUM_ID}.recal_1.failed ; 
aws s3 cp ${TUM_ID}.recal_data.grp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/4_split_bams/ ; 
rm ${TUM_ID}.split.bam; 
rm ${TUM_ID}.split.bai; 
rm ${TUM_ID}.recal_data.grp;
echo `date "+%Y-%m-%d %H:%M:%S"` "BaseRecalibrator done" >> log.txt

## base quality recalibration - apply BQSR
echo `date "+%Y-%m-%d %H:%M:%S"` "ApplyBQSR start" >> log.txt
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/4_split_bams/${TUM_ID}.split.bam ./ ; 
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/4_split_bams/${TUM_ID}.split.bai ./ ;
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/4_split_bams/${TUM_ID}.recal_data.grp ./ ; 
/home/ubuntu/tools/gatk-4.4.0.0/gatk ApplyBQSR \
-I ${TUM_ID}.split.bam \
-O ${TUM_ID}.recal.bam \
-R /stereoseq/reference/GRCh38.109.fa \
-bqsr ${TUM_ID}.recal_data.grp \
--tmp-dir /stereoseq/tmp \
&& touch ${TUM_ID}.recal_2.success || touch ${TUM_ID}.recal_2.failed ; 
aws s3 cp ${TUM_ID}.recal.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/5_recal_bams/ ; 
aws s3 cp ${TUM_ID}.recal.bai s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/bams/5_recal_bams/ ; 
rm ${TUM_ID}.split.bam; 
rm ${TUM_ID}.split.bai; 
rm ${TUM_ID}.recal.bam; 
rm ${TUM_ID}.recal.bai; 
rm ${TUM_ID}.recal_data.grp;
echo `date "+%Y-%m-%d %H:%M:%S"` "ApplyBQSR done" >> log.txt