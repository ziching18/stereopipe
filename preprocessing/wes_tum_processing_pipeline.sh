# TUMOUR PRE-PROCESSING (FQ) - WES
# while read line; do slx=$(echo $line | awk '{print $1}'); id=$(echo $line | awk '{print $2}'); cd /stereoseq/all_samples/normal/$id; bash /stereoseq/code/fastq_preprocessing_pipeline.sh $slx $id; done < /stereoseq/code/samples.txt

# # path slx slx_id id
# crm.tumorstudy.ngs.raw/	14336	14336.i710_i504	SD0693
# crm.tumorstudy.ngs.raw/	14336	14336.i709_i505	SD0507
# crm.tumorstudy.ngs.raw/	14499	14499.i704_i503	SD1043
# crm.tumorstudy.ngs.raw/	14499	14499.i702_i507	SD0560
# crm.tumorstudy.ngs.raw/	15093	15093.i704_i502	SD0683
# crm.tumorstudy.ngs.raw/WES/	15093	15093.i706_i517	SD1182
# crm.tumorstudy.ngs.raw/WES/	15098	15098.i707_i517	SD0781
# crm.tumorstudy.ngs.raw/WES/	15098	15098.i710_i517	SD1225
# crm.tumorstudy.ngs.raw/WES/	15098	15098.i707_i517	SD1043

CPATH=$1;
SLX=$2; 
SLX_ID=$3; 
TUM_ID=$4; 
cd /stereoseq/all_samples/wes/${TUM_ID}/;
### trim
#ls |
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} Trim Galore start" >> log.txt
##SLX-14290.i708_i508.HL2T2BBXX.s_8.r_1.fq.gz
aws s3 ls ${CPATH}/SLX-${SLX}/ | grep SLX-${SLX_ID}. |
awk '{print $NF}'|
grep 'fq.gz$' | grep -v 0000 |
cut -d '.' -f 1,2,3,4 | 
sort | uniq | 
xargs -P 2 -I AAA sh -c  \
" aws s3 cp s3://${CPATH}/SLX-${SLX}/AAA.r_1.fq.gz ./;
aws s3 cp s3://${CPATH}/SLX-${SLX}/AAA.r_2.fq.gz ./;
trim_galore AAA.r_1.fq.gz AAA.r_2.fq.gz --paired --gzip -o . -a CTGTCTCTTATACACATCT 2>AAA.log ; \
aws s3 cp AAA.r_1_val_1.fq.gz s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/1_trim_galore_out/ ; \
aws s3 cp AAA.r_2_val_2.fq.gz s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/1_trim_galore_out/ && touch AAA.trim.success || touch AAA.trim.failed ; \
rm AAA.*.fq.gz ; \
rm AAA.*.bam ; \
rm AAA.*.fq " ;
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} Trim Galore finish" >> log.txt

### map
#ls |
cd /stereoseq/all_samples/wes/${TUM_ID}/;
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} BWA start" >> log.txt
file=$(echo $SLX_ID | cut -d '.' -f 1)
#echo ${file} | xargs -P 4 -I %% sh -c "aws s3 ls crm.sequencing.raw.data.sharing/batch1/SLX-%%/" | grep SLX-${SLX_ID} |
aws s3 ls s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/1_trim_galore_out/ | #grep SLX-${SLX_ID}. |
awk '{print $NF}' |
grep 'fq.gz$' | grep -v 0000 |
cut -d '.' -f 1,2,3,4 | 
sort | uniq | 
xargs -P 2 -I AAA sh -c  \
" aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/1_trim_galore_out/AAA.r_1_val_1.fq.gz ./ ; \
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/1_trim_galore_out/AAA.r_2_val_2.fq.gz ./ ; \
zcat AAA.r_1_val_1.fq.gz | awk '(NR%2==0){\$0=substr(\$0,1,75)}{print}' > AAA.r_1_bwa_in.fq ; \
zcat AAA.r_2_val_2.fq.gz | awk '(NR%2==0){\$0=substr(\$0,1,75)}{print}' > AAA.r_2_bwa_in.fq ; \
bwa mem -M -t 4 /stereoseq/reference/GRCh38.109.bwa.fa \
AAA.r_1_bwa_in.fq \
AAA.r_2_bwa_in.fq 2>AAA.bwa.log | \
samtools view --threads 8 -b - | \
samtools sort --threads 8 > AAA.sorted.bam ; \
aws s3 cp AAA.bwa.log s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/2_bwa_out/ ; \
aws s3 cp AAA.sorted.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/2_bwa_out/ && touch AAA.map.success || touch AAA.map.failed ; \
rm AAA.*.fq.gz; \
rm AAA.*.fq; " ;
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} BWA finish" >> log.txt

## merge
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} MergeSamFiles start" >> log.txt
cd /stereoseq/all_samples/wes/${TUM_ID}/;
aws s3 ls crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/2_bwa_out/ | grep SLX-${SLX_ID}. |
grep .sorted.bam$ | awk '{print $NF}'| while read l; do
    echo "-I "$l
done > ${TUM_ID}_wes_tum_bams.list
/home/ubuntu/tools/gatk-4.4.0.0/gatk MergeSamFiles --USE_THREADING true \
--arguments_file ${TUM_ID}_wes_tum_bams.list \
-O ${TUM_ID}.wes_tum.merged.bam
samtools index ${TUM_ID}.wes_tum.merged.bam ${TUM_ID}.wes_tum.merged.bai;
aws s3 cp ${TUM_ID}.wes_tum.merged.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/3_merged_bams/; 
aws s3 cp ${TUM_ID}.wes_tum.merged.bai s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/3_merged_bams/;
rm *.bam;
rm *.bai;
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} MergeSamFiles finish" >> log.txt

### add read groups
#ls |
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} AddOrReplaceReadGroups start" >> log.txt
cd /stereoseq/all_samples/wes/${TUM_ID}/;
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/3_merged_bams/${TUM_ID}.wes_tum.merged.bam ./; 
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/3_merged_bams/${TUM_ID}.wes_tum.merged.bai ./; 
/home/ubuntu/tools/gatk-4.4.0.0/gatk AddOrReplaceReadGroups \
-I ${TUM_ID}.wes_tum.merged.bam \
-O ${TUM_ID}.wes_tum.merged.RGA.bam \
--RGID ${SLX_ID} \
--RGLB $(echo ${SLX_ID} | cut -d '.' -f 1) \
--RGPL Illumina \
--RGPU ${SLX_ID}.${TUM_ID} \
--RGSM ${TUM_ID}_MN \
&& touch ${TUM_ID}.addr.success || touch ${TUM_ID}.addr.failed ; 
aws s3 cp ${TUM_ID}.wes_tum.merged.RGA.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/4_rga_out/ ; 
rm ${TUM_ID}.wes_tum.merged.bam; 
rm ${TUM_ID}.wes_tum.merged.bai;
rm ${TUM_ID}.wes_tum.merged.RGA.bam ; 
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} AddOrReplaceReadGroups finish" >> log.txt

## mark duplicates
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} MarkDuplicates start" >> log.txt
cd /stereoseq/all_samples/wes/${TUM_ID}/;
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/4_rga_out/${TUM_ID}.wes_tum.merged.RGA.bam ./ ; 
/home/ubuntu/tools/gatk-4.4.0.0/gatk MarkDuplicates \
-I ${TUM_ID}.wes_tum.merged.RGA.bam \
-O ${TUM_ID}.wes_tum.dedup.bam \
-M ${TUM_ID}.wes_tum.MarkDup.metrics \
--CREATE_INDEX true \
--VALIDATION_STRINGENCY SILENT ; 
aws s3 cp ${TUM_ID}.wes_tum.dedup.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/5_dedupped_bams/ ; 
aws s3 cp ${TUM_ID}.wes_tum.dedup.bai s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/5_dedupped_bams/ ; 
aws s3 cp ${TUM_ID}.wes_tum.MarkDup.metrics s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/5_dedupped_bams/ ; 
rm ${TUM_ID}.wes_tum.dedup.bam; 
rm ${TUM_ID}.wes_tum.dedup.bai; 
rm ${TUM_ID}.wes_tum.merged.RGA.bam; 
rm ${TUM_ID}.wes_tum.MarkDup.metrics ;
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} MarkDuplicates finish" >> log.txt


## split reads
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} SplitNCigarReads start" >> log.txt
cd /stereoseq/all_samples/wes/${TUM_ID}/;
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/5_dedupped_bams/${TUM_ID}.wes_tum.dedup.bam ./ ; 
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/5_dedupped_bams/${TUM_ID}.wes_tum.dedup.bai ./ ;  
/home/ubuntu/tools/gatk-4.4.0.0/gatk SplitNCigarReads \
-I ${TUM_ID}.wes_tum.dedup.bam \
-O ${TUM_ID}.wes_tum.split.bam \
-R /stereoseq/reference/GRCh38.109.fa \
--tmp-dir /stereoseq/tmp \
&& touch ${TUM_ID}.wes_tum.split.success || touch ${TUM_ID}.wes_tum.split.failed ; 
aws s3 cp ${TUM_ID}.wes_tum.split.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/6_split_bams/ ; 
aws s3 cp ${TUM_ID}.wes_tum.split.bai s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/6_split_bams/ ; 
rm ${TUM_ID}.wes_tum.dedup.bam; 
rm ${TUM_ID}.wes_tum.dedup.bai; 
rm ${TUM_ID}.wes_tum.split.bam; 
rm ${TUM_ID}.wes_tum.split.bai; 
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} SplitNCigarReads finish" >> log.txt

## base quality recalibration - split
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} BaseRecalibrator start" >> log.txt
cd /stereoseq/all_samples/wes/${TUM_ID}/;
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/6_split_bams/${TUM_ID}.wes_tum.split.bam ./ ; 
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/6_split_bams/${TUM_ID}.wes_tum.split.bai ./ ;  
/home/ubuntu/tools/gatk-4.4.0.0/gatk BaseRecalibrator \
-I ${TUM_ID}.wes_tum.split.bam \
-O ${TUM_ID}.wes_tum.recal_data.grp \
-R /stereoseq/reference/GRCh38.109.fa \
--tmp-dir /stereoseq/tmp \
--known-sites /stereoseq/reference/Homo_sapiens_assembly38.known_indels.renamed.vcf \
--known-sites /stereoseq/reference/Homo_sapiens_assembly38.dbsnp138.renamed.vcf \
--known-sites /stereoseq/reference/Mills_and_1000G_gold_standard.indels.hg38.renamed.vcf \
&& touch ${TUM_ID}.wes_tum.recal_1.success || touch ${TUM_ID}.wes_tum.recal_1.failed ; 
aws s3 cp ${TUM_ID}.wes_tum.recal_data.grp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/6_split_bams/ ; 
rm ${TUM_ID}.wes_tum.split.bam; 
rm ${TUM_ID}.wes_tum.split.bai; 
rm ${TUM_ID}.wes_tum.recal_data.grp ; 
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} BaseRecalibrator finish" >> log.txt

## base quality recalibration - apply BQSR
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} ApplyBQSR start" >> log.txt
cd /stereoseq/all_samples/wes/${TUM_ID}/;
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/6_split_bams/${TUM_ID}.wes_tum.split.bam ./ ;  
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/6_split_bams/${TUM_ID}.wes_tum.split.bai ./ ; 
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/6_split_bams/${TUM_ID}.wes_tum.recal_data.grp ./ ; 
/home/ubuntu/tools/gatk-4.4.0.0/gatk ApplyBQSR \
-I ${TUM_ID}.wes_tum.split.bam \
-O ${TUM_ID}.wes_tum.recal.bam \
-R /stereoseq/reference/GRCh38.109.fa \
-bqsr ${TUM_ID}.wes_tum.recal_data.grp \
--tmp-dir /stereoseq/tmp \
&& touch ${TUM_ID}.wes_tum.recal_2.success || touch ${TUM_ID}.wes_tum.recal_2.failed ; 
aws s3 cp ${TUM_ID}.wes_tum.recal.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/7_recal_bams/ ; 
aws s3 cp ${TUM_ID}.wes_tum.recal.bai s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum/7_recal_bams/ ; 
rm ${TUM_ID}.wes_tum.split.bam; 
rm ${TUM_ID}.wes_tum.split.bai; 
rm ${TUM_ID}.wes_tum.recal.bam; 
rm ${TUM_ID}.wes_tum.recal.bai; 
rm ${TUM_ID}.wes_tum.recal_data.grp; 
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} ApplyBQSR finish" >> log.txt

