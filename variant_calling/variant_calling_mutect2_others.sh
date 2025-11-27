## VARIANT CALLING W/ MUTECT2

#SLX_ID=$1; 
TUM_ID=$1; 
echo ${TUM_ID}
## calling
#aws s3 sync s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal/6_recal_bams/ /stereoseq/all_samples/${TUM_ID}/normal/;

#cd /stereoseq/all_samples/normal/${TUM_ID}/;

# make name list of normal
#aws s3 ls s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/normal_updatedFeb25/7_recal_bams/ | grep SLX-${SLX_ID}. |
#grep .bam$ | awk '{print $NF}' | cut -d '.' -f 1,2,3,4 | while read l; do
#    echo "-normal "$l
#done > ${TUM_ID}_normal_names.list

cd /stereoseq/all_samples/wes/${TUM_ID}/;

aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum_corrected/7_recal_bams/${TUM_ID}.wes_tum.recal.bam ./;
aws s3 cp s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum_corrected/7_recal_bams/${TUM_ID}.wes_tum.recal.bai ./;

/home/ubuntu/tools/gatk-4.4.0.0/gatk Mutect2 \
-R /stereoseq/reference/GRCh38.109.fa \
-I /stereoseq/all_samples/wes/${TUM_ID}/${TUM_ID}.wes_tum.recal.bam \
-I /stereoseq/all_samples/normal/${TUM_ID}/${TUM_ID}.normal.recal.bam \
--panel-of-normals /stereoseq/all_samples/normal/${TUM_ID}/${TUM_ID}_pon.vcf.gz \
-normal ${TUM_ID}_MN \
--tmp-dir /stereoseq/tmp \
-O ${TUM_ID}.wes_tum.somatic.vcf.gz \
2>${TUM_ID}.wes_tum.somatic.vcf.gz.log ; 
gzip -dkf ${TUM_ID}.wes_tum.somatic.vcf.gz; # unzip
aws s3 cp ${TUM_ID}.wes_tum.somatic.vcf.gz s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum_corrected/8_vcf/;
aws s3 cp ${TUM_ID}.wes_tum.somatic.vcf.gz.log s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum_corrected/8_vcf/;
aws s3 cp ${TUM_ID}.wes_tum.somatic.vcf s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/wes_tum_corrected/8_vcf/;
#--arguments_file /stereoseq/all_samples/normal/${TUM_ID}/${TUM_ID}_normal_names.list \

rm ${TUM_ID}.wes_tum.recal.bam;
rm ${TUM_ID}.wes_tum.recal.bai;

## variant filtration pt. 1
/home/ubuntu/tools/gatk-4.4.0.0/gatk FilterMutectCalls \
--variant ${TUM_ID}.wes_tum.somatic.vcf.gz \
--output ${TUM_ID}.wes_tum.somatic.mutect2filtered.vcf.gz \
--reference /stereoseq/reference/GRCh38.109.fa
gzip -dkf ${TUM_ID}.wes_tum.somatic.mutect2filtered.vcf.gz

### see filtration results, we want PASS
grep -v "^#" ${TUM_ID}.wes_tum.somatic.mutect2filtered.vcf | cut -f7 | sort | uniq -c > ${TUM_ID}.wes_tum.filtered_calls.txt
### filter for PASS
bcftools view -f PASS ${TUM_ID}.wes_tum.somatic.mutect2filtered.vcf -Oz -o ${TUM_ID}.wes_tum.somatic.filtered1.vcf 

## funcotator
/home/ubuntu/tools/gatk-4.4.0.0/gatk Funcotator \
--variant ${TUM_ID}.wes_tum.somatic.filtered1.vcf \
--reference /stereoseq/reference/GRCh38.109.fa \
--ref-version hg38 \
--data-sources-path /stereoseq/reference/funcotator/funcotator_dataSources.v1.7.20200521s \
--output ${TUM_ID}.wes_tum.somatic.filtered1.funcotated.vcf \
--output-file-format VCF;
#aws s3 cp ${TUM_ID}.somatic.filtered.funcotated.vcf s3://crm.steroseq.raw.data/Breast_CACRMY/all_samples/${TUM_ID}/vcf_updatedFeb25/;

## variant filtration pt. 2
### VAF < 0.1? + < 3 supporting reads
### VAF > 0.8
bcftools view -e 'FORMAT/AF[1:0] > 0.8' ${TUM_ID}.wes_tum.somatic.filtered1.funcotated.vcf -Oz -o ${TUM_ID}.wes_tum.somatic.filtered2.funcotated.vcf
### TLOD < 5.6
bcftools view -e 'INFO/TLOD < 5.6' ${TUM_ID}.wes_tum.somatic.filtered2.funcotated.vcf -Oz -o ${TUM_ID}.wes_tum.somatic.filtered3.funcotated.vcf
### no short tandem repeats (STR)
bcftools view -e 'INFO/STR=1' ${TUM_ID}.wes_tum.somatic.filtered3.funcotated.vcf -Oz -o ${TUM_ID}.wes_tum.somatic.filtered4.funcotated.vcf

## variant filtration pt. 3
grep -v "#" ${TUM_ID}.wes_tum.somatic.filtered4.funcotated.vcf | \
awk '$1 <= 23 && $1 >=1' > ${TUM_ID}.wes_tum.somatic.filtered5.funcotated.vcf
# #grep -v "\[HLA" \
grep "#" ${TUM_ID}.wes_tum.somatic.filtered4.funcotated.vcf > ${TUM_ID}.wes_tum.somatic.filtered6.funcotated.vcf
awk '$1 <= 23 && $1 >=1' ${TUM_ID}.wes_tum.somatic.filtered4.funcotated.vcf >> ${TUM_ID}.wes_tum.somatic.filtered6.funcotated.vcf

## variant filtration pt. 4
grep -v -e "\[IGH" -e "\[IGK" -e "\[IGL" -e "\[HLA" ${TUM_ID}.wes_tum.somatic.filtered6.funcotated.vcf \
> ${TUM_ID}.wes_tum.somatic.filtered7.funcotated.vcf
#| awk '{print $8}' | cut -d ';' -f 5 |  cut -d '|' -f 1 | cut -d '[' -f 2 > IG.txt
