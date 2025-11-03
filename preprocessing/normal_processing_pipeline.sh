# MATCHED NORMAL PRE-PROCESSING (FQ)
# while read line; do slx=$(echo $line | awk '{print $1}'); id=$(echo $line | awk '{print $2}'); cd /stereoseq/all_samples/normal/$id; bash /stereoseq/code/fastq_preprocessing_pipeline.sh $slx $id; done < /stereoseq/code/samples.txt


SLX_ID=$1; 
TUM_ID=$2; 
cd /stereoseq/all_samples/normal/${TUM_ID}/;
### trim
#ls |
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} Trim Galore start" >> log.txt
file=$(echo $SLX_ID | cut -d '.' -f 1)
echo ${file} | xargs -P 4 -I %% sh -c "aws s3 ls crm.sequencing.raw.data.sharing/batch1/SLX-%%/" | grep SLX-${SLX_ID}. |
awk '{print $NF}'|
grep 'fq.gz$' | grep -v 0000 |
cut -d '.' -f 1,2,3,4 | 
sort | uniq | 
xargs -P 2 -I AAA sh -c  \
" aws s3 cp s3://crm.sequencing.raw.data.sharing/batch1/SLX-${file}/AAA.r_1.fq.gz ./;
aws s3 cp s3://crm.sequencing.raw.data.sharing/batch1/SLX-${file}/AAA.r_2.fq.gz ./;
trim_galore AAA.r_1.fq.gz AAA.r_2.fq.gz --paired --gzip -o . -a CTGTCTCTTATACACATCT 2>AAA.log ; \
aws s3 cp AAA.r_1_val_1.fq.gz s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/1_trim_galore_out/ ; \
aws s3 cp AAA.r_2_val_2.fq.gz s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/1_trim_galore_out/ && touch AAA.trim.success || touch AAA.trim.failed ; \
rm AAA.*.fq.gz ; \
rm AAA.*.bam ; \
rm AAA.*.fq " ;
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} Trim Galore finish" >> log.txt

### map
#ls |
cd /stereoseq/all_samples/normal/${TUM_ID}/;
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} BWA start" >> log.txt
file=$(echo $SLX_ID | cut -d '.' -f 1)
#echo ${file} | xargs -P 4 -I %% sh -c "aws s3 ls crm.sequencing.raw.data.sharing/batch1/SLX-%%/" | grep SLX-${SLX_ID} |
aws s3 ls s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/1_trim_galore_out/ | #grep SLX-${SLX_ID}. |
awk '{print $NF}' |
grep 'fq.gz$' | grep -v 0000 |
cut -d '.' -f 1,2,3,4 | 
sort | uniq | 
xargs -P 2 -I AAA sh -c  \
" aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/1_trim_galore_out/AAA.r_1_val_1.fq.gz ./ ; \
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/1_trim_galore_out/AAA.r_2_val_2.fq.gz ./ ; \
zcat AAA.r_1_val_1.fq.gz | awk '(NR%2==0){\$0=substr(\$0,1,75)}{print}' > AAA.r_1_bwa_in.fq ; \
zcat AAA.r_2_val_2.fq.gz | awk '(NR%2==0){\$0=substr(\$0,1,75)}{print}' > AAA.r_2_bwa_in.fq ; \
bwa mem -M -t 4 /stereoseq/reference/GRCh38.109.bwa.fa \
AAA.r_1_bwa_in.fq \
AAA.r_2_bwa_in.fq 2>AAA.bwa.log | \
samtools view --threads 8 -b - | \
samtools sort --threads 8 > AAA.sorted.bam ; \
aws s3 cp AAA.bwa.log s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/2_bwa_out/ ; \
aws s3 cp AAA.sorted.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/2_bwa_out/ && touch AAA.map.success || touch AAA.map.failed ; \
rm AAA.*.bam; \
rm AAA.*.fq.gz; \
rm AAA.*.fq; " ;
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} BWA finish" >> log.txt

## merge
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} MergeSamFiles start" >> log.txt
cd /stereoseq/all_samples/normal/${TUM_ID}/;
aws s3 ls crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/2_bwa_out/ | grep SLX-${SLX_ID}. |
grep .sorted.bam$ | awk '{print $NF}'| while read l; do
    echo "-I "$l
done > ${TUM_ID}_normal_bams.list
/home/ubuntu/tools/gatk-4.4.0.0/gatk MergeSamFiles --USE_THREADING true \
--arguments_file ${TUM_ID}_normal_bams.list \
-O ${TUM_ID}.normal.merged.bam
samtools index ${TUM_ID}.normal.merged.bam ${TUM_ID}.normal.merged.bai;
aws s3 cp ${TUM_ID}.normal.merged.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/3_merged_bams/; 
aws s3 cp ${TUM_ID}.normal.merged.bai s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/3_merged_bams/;
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} MergeSamFiles finish" >> log.txt

### add read groups
#ls |
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} AddOrReplaceReadGroups start" >> log.txt
cd /stereoseq/all_samples/normal/${TUM_ID}/;
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/3_merged_bams/${TUM_ID}.normal.merged.bam ./; 
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/3_merged_bams/${TUM_ID}.normal.merged.bai ./; 
/home/ubuntu/tools/gatk-4.4.0.0/gatk AddOrReplaceReadGroups \
-I ${TUM_ID}.normal.merged.bam \
-O ${TUM_ID}.normal.merged.RGA.bam \
--RGID ${SLX_ID} \
--RGLB $(echo ${SLX_ID} | cut -d '.' -f 1) \
--RGPL Illumina \
--RGPU ${SLX_ID}.${TUM_ID} \
--RGSM ${TUM_ID}_MN \
&& touch ${TUM_ID}.addr.success || touch ${TUM_ID}.addr.failed ; 
aws s3 cp ${TUM_ID}.normal.merged.RGA.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/4_rga_out/ ; 
rm ${TUM_ID}.normal.merged.bam; 
rm ${TUM_ID}.normal.merged.bai;
rm ${TUM_ID}.normal.merged.RGA.bam ; 
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} AddOrReplaceReadGroups finish" >> log.txt

## mark duplicates
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} MarkDuplicates start" >> log.txt
cd /stereoseq/all_samples/normal/${TUM_ID}/;
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/4_rga_out/${TUM_ID}.normal.merged.RGA.bam ./ ; 
/home/ubuntu/tools/gatk-4.4.0.0/gatk MarkDuplicates \
-I ${TUM_ID}.normal.merged.RGA.bam \
-O ${TUM_ID}.normal.dedup.bam \
-M ${TUM_ID}.normal.MarkDup.metrics \
--CREATE_INDEX true \
--VALIDATION_STRINGENCY SILENT ; 
aws s3 cp ${TUM_ID}.normal.dedup.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/5_dedupped_bams/ ; 
aws s3 cp ${TUM_ID}.normal.dedup.bai s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/5_dedupped_bams/ ; 
aws s3 cp ${TUM_ID}.normal.MarkDup.metrics s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/5_dedupped_bams/ ; 
rm ${TUM_ID}.normal.dedup.bam; 
rm ${TUM_ID}.normal.dedup.bai; 
rm ${TUM_ID}.normal.merged.RGA.bam; 
rm ${TUM_ID}.normal.MarkDup.metrics ;
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} MarkDuplicates finish" >> log.txt


## split reads
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} SplitNCigarReads start" >> log.txt
cd /stereoseq/all_samples/normal/${TUM_ID}/;
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/5_dedupped_bams/${TUM_ID}.normal.dedup.bam ./ ; 
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/5_dedupped_bams/${TUM_ID}.normal.dedup.bai ./ ;  
/home/ubuntu/tools/gatk-4.4.0.0/gatk SplitNCigarReads \
-I ${TUM_ID}.normal.dedup.bam \
-O ${TUM_ID}.normal.split.bam \
-R /stereoseq/reference/GRCh38.109.fa \
--tmp-dir /stereoseq/tmp \
&& touch ${TUM_ID}.normal.split.success || touch ${TUM_ID}.normal.split.failed ; 
aws s3 cp ${TUM_ID}.normal.split.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/6_split_bams/ ; 
aws s3 cp ${TUM_ID}.normal.split.bai s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/6_split_bams/ ; 
rm ${TUM_ID}.normal.dedup.bam; 
rm ${TUM_ID}.normal.dedup.bai; 
rm ${TUM_ID}.normal.split.bam; 
rm ${TUM_ID}.normal.split.bai; 
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} SplitNCigarReads finish" >> log.txt

## base quality recalibration - split
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} BaseRecalibrator start" >> log.txt
cd /stereoseq/all_samples/normal/${TUM_ID}/;
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/6_split_bams/${TUM_ID}.normal.split.bam ./ ; 
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/6_split_bams/${TUM_ID}.normal.split.bai ./ ;  
/home/ubuntu/tools/gatk-4.4.0.0/gatk BaseRecalibrator \
-I ${TUM_ID}.normal.split.bam \
-O ${TUM_ID}.normal.recal_data.grp \
-R /stereoseq/reference/GRCh38.109.fa \
--tmp-dir /stereoseq/tmp \
--known-sites /stereoseq/reference/Homo_sapiens_assembly38.known_indels.renamed.vcf \
--known-sites /stereoseq/reference/Homo_sapiens_assembly38.dbsnp138.renamed.vcf \
--known-sites /stereoseq/reference/Mills_and_1000G_gold_standard.indels.hg38.renamed.vcf \
&& touch ${TUM_ID}.normal.recal_1.success || touch ${TUM_ID}.normal.recal_1.failed ; 
aws s3 cp ${TUM_ID}.normal.recal_data.grp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/6_split_bams/ ; 
rm ${TUM_ID}.normal.split.bam; 
rm ${TUM_ID}.normal.split.bai; 
rm ${TUM_ID}.normal.recal_data.grp ; 
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} BaseRecalibrator finish" >> log.txt

## base quality recalibration - apply BQSR
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} ApplyBQSR start" >> log.txt
cd /stereoseq/all_samples/normal/${TUM_ID}/;
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/6_split_bams/${TUM_ID}.normal.split.bam ./ ;  
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/6_split_bams/${TUM_ID}.normal.split.bai ./ ; 
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/6_split_bams/${TUM_ID}.normal.recal_data.grp ./ ; 
/home/ubuntu/tools/gatk-4.4.0.0/gatk ApplyBQSR \
-I ${TUM_ID}.normal.split.bam \
-O ${TUM_ID}.normal.recal.bam \
-R /stereoseq/reference/GRCh38.109.fa \
-bqsr ${TUM_ID}.normal.recal_data.grp \
--tmp-dir /stereoseq/tmp \
&& touch ${TUM_ID}.normal.recal_2.success || touch ${TUM_ID}.normal.recal_2.failed ; 
aws s3 cp ${TUM_ID}.normal.recal.bam s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/7_recal_bams/ ; 
aws s3 cp ${TUM_ID}.normal.recal.bai s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/7_recal_bams/ ; 
rm ${TUM_ID}.normal.split.bam; 
rm ${TUM_ID}.normal.split.bai; 
rm ${TUM_ID}.normal.recal.bam; 
rm ${TUM_ID}.normal.recal.bai; 
rm ${TUM_ID}.normal.recal_data.grp; 
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} ApplyBQSR finish" >> log.txt

## genotyping w/ haplotypecaller - calling all variants
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} HaplotypeCaller start" >> log.txt
cd /stereoseq/all_samples/normal/${TUM_ID}/;
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/7_recal_bams/${TUM_ID}.normal.recal.bam ./ ;  
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/7_recal_bams/${TUM_ID}.normal.recal.bai ./ ;
/home/ubuntu/tools/gatk-4.4.0.0/gatk HaplotypeCaller \
-R /stereoseq/reference/GRCh38.109.fa \
-I /stereoseq/all_samples/normal/${TUM_ID}/${TUM_ID}.normal.recal.bam \
-O ${TUM_ID}.normal.all.vcf.gz \
--tmp-dir /stereoseq/tmp \
--dont-use-soft-clipped-bases \
2>${TUM_ID}.normal.all.vcf.gz.log ; 
gzip -dkf ${TUM_ID}.normal.all.vcf.gz; # unzip
aws s3 cp ${TUM_ID}.normal.all.vcf.gz s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/vcf_updatedFeb25/;
aws s3 cp ${TUM_ID}.normal.all.vcf.gz.log s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/vcf_updatedFeb25/;
aws s3 cp ${TUM_ID}.normal.all.vcf s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/vcf_updatedFeb25/;
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} HaplotypeCaller finish" >> log.txt

### Mutect2 ###
#ls |
# make panel of normals (1)
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} Mutect2 make panel of normals start" >> log.txt
cd /stereoseq/all_samples/normal/${TUM_ID}/;
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/7_recal_bams/${TUM_ID}.normal.recal.bam ./ ; 
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/7_recal_bams/${TUM_ID}.normal.recal.bai ./ ; 
/home/ubuntu/tools/gatk-4.4.0.0/gatk Mutect2 \
-I ${TUM_ID}.normal.recal.bam \
-O ${TUM_ID}.normal.vcf.gz \
-R /stereoseq/reference/GRCh38.109.fa \
--max-mnp-distance 0 \
2>${TUM_ID}.normal.vcf.gz.log; 
aws s3 cp ${TUM_ID}.normal.vcf.gz s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/8_mutect2_pon/ ; 
aws s3 cp ${TUM_ID}.normal.vcf.gz.tbi s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/8_mutect2_pon/ ; 
rm ${TUM_ID}.normal.recal.bam; 
rm ${TUM_ID}.normal.recal.bai ; 
echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} Mutect2 make panel of normals finish" >> log.txt


# make panel of normals (2)
# gatk Mutect2 -R reference.fasta -I normal1.bam --max-mnp-distance 0 -O normal1.vcf.gz
## combine into one big matched normal

aws s3 ls crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/8_mutect2_pon/ | 
grep vcf.gz$ | awk '{print $NF}'| while read l; do
    echo "-V "$l
done > ${TUM_ID}_PON_arg.list

/home/ubuntu/tools/gatk-4.4.0.0/gatk GenomicsDBImport \
-R /stereoseq/reference/GRCh38.109.fa \
-L /stereoseq/reference/wgs_calling_regions.hg38.renamed.interval_list \
--genomicsdb-workspace-path ${TUM_ID}_pon_db \
--tmp-dir /stereoseq/tmp \
-V ${TUM_ID}.normal.vcf.gz ;
#--arguments_file ${TUM_ID}_PON_arg.list ;
#-L /stereoseq/reference/hg38.even.handcurated.20k.intervals \

/home/ubuntu/tools/gatk-4.4.0.0/gatk CreateSomaticPanelOfNormals \
-R /stereoseq/reference/GRCh38.109.fa \
--germline-resource /stereoseq/reference/Mills_and_1000G_gold_standard.indels.hg38.renamed.vcf \
-V gendb://${TUM_ID}_pon_db \
-O ${TUM_ID}_pon.vcf.gz \
2>${TUM_ID}_normal_2ormo_PON.vcf.gz.log;
aws s3 cp ${TUM_ID}_pon.vcf.gz s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/8_mutect2_pon/ ;
aws s3 cp ${TUM_ID}_normal_2ormo_PON.vcf.gz.log s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/8_mutect2_pon/ ;


