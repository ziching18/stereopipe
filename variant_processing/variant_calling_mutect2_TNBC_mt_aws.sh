## MITOCHONDRIAL VARIANT CALLING W/ MUTECT2 ON WES BAMS IN AWS
n=$1;

# parallel bash /stereoseq/code/stereopipe/variant_processing/variant_calling_mutect2_TNBC_mt_aws.sh ::: $(seq 6)

grep Tumor_Plate${n} /stereoseq/wes/MyBrCa_WES_TNBC_ids.txt | 
awk '{print $1}' |
xargs -I % sh -c \
"

aws s3 cp s3://crm.tumorstudy.mamduh/WES/Tumor_Plate${n}/BQSR/%_TUM.recalibrated.bam /stereoseq/wes/bam/plate${n}/ --force-glacier-transfer; \
aws s3 cp s3://crm.tumorstudy.mamduh/WES/Tumor_Plate${n}/BQSR/%_TUM.recalibrated.bai /stereoseq/wes/bam/plate${n}/ --force-glacier-transfer; \


samtools view -b /stereoseq/wes/bam/plate${n}/%.recalibrated.bam MT > /stereoseq/wes/bam/plate${n}/%.recalibrated.MT.bam; \
samtools index /stereoseq/wes/bam/plate${n}/%.recalibrated.MT.bam; \


~/tools/gatk-4.6.2.0/gatk Mutect2 \
   -R /stereoseq/reference/hs37d5.fa.gz \
   -I /stereoseq/wes/bam/plate${n}/%.recalibrated.MT.bam \
   --mitochondria-mode true \
   -O /stereoseq/wes/vcf/plate${n}/%.MT.vcf; \


~/tools/gatk-4.6.2.0/gatk FilterMutectCalls \
   -V /stereoseq/wes/vcf/plate${n}/%.MT.vcf \
   -R /stereoseq/reference/hs37d5.fa.gz \
   -O /stereoseq/wes/vcf/plate${n}/%.MT.filtered.vcf; \


bcftools view -f PASS /stereoseq/wes/vcf/plate${n}/%.MT.filtered.vcf -Oz -o /stereoseq/wes/vcf/plate${n}/%.MT.filtered.PASS.vcf; \


aws s3 cp /stereoseq/wes/vcf/plate${n}/%.MT.filtered.vcf s3://crm.tumorstudy.mamduh/WES/somatic_variants/GATK4_Mutect2_MT_variants/%.MT.filtered.vcf; \
aws s3 cp /stereoseq/wes/vcf/plate${n}/%.MT.filtered.vcf.idx s3://crm.tumorstudy.mamduh/WES/somatic_variants/GATK4_Mutect2_MT_variants/%.MT.filtered.vcf.idx; \
aws s3 cp /stereoseq/wes/vcf/plate${n}/%.MT.filtered.PASS.vcf s3://crm.tumorstudy.mamduh/WES/somatic_variants/GATK4_Mutect2_MT_variants/%.MT.filtered.PASS.vcf; \
aws s3 cp /stereoseq/wes/vcf/plate${n}/%.MT.filtered.PASS.vcf.idx s3://crm.tumorstudy.mamduh/WES/somatic_variants/GATK4_Mutect2_MT_variants/%.MT.filtered.PASS.vcf.idx; \


rm /stereoseq/wes/bam/plate${n}/*;
rm /stereoseq/wes/vcf/plate${n}/*;

";
