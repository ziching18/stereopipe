## MITOCHONDRIAL VARIANT CALLING W/ MUTECT2 ON WES BAMS IN AWS

n=$1

aws s3 ls s3://crm.tumorstudy.mamduh/WES/Tumor_Plate${n}/BQSR/ |
awk '{print $NF}'|
grep 'bam' | 
cut -d '.' -f 1 | 
sort | uniq | 
xargs -I % sh -c \
"

aws s3 cp s3://crm.tumorstudy.mamduh/WES/Tumor_Plate${n}/BQSR/%.recalibrated.bam /stereoseq/wes/bam/; \
aws s3 cp s3://crm.tumorstudy.mamduh/WES/Tumor_Plate${n}/BQSR/%.recalibrated.bai /stereoseq/wes/bam/; \


samtools view -b /stereoseq/wes/bam/%.recalibrated.bam MT > /stereoseq/wes/bam/%.recalibrated.MT.bam; \
samtools index /stereoseq/wes/bam/%.recalibrated.MT.bam; \


aws s3 cp /stereoseq/wes/bam/%.recalibrated.MT.bam s3://crm.tumorstudy.mamduh/WES/MT/%.recalibrated.MT.bam; \


~/tools/gatk-4.6.2.0/gatk Mutect2 \
   -R /stereoseq/reference/hs37d5.fa.gz \
   -I /stereoseq/wes/bam/%.recalibrated.MT.bam \
   --mitochondria-mode true \
   -O /stereoseq/wes/vcf/%.MT.vcf; \


~/tools/gatk-4.6.2.0/gatk FilterMutectCalls \
   -V /stereoseq/wes/vcf/%.MT.vcf \
   -R /stereoseq/reference/hs37d5.fa.gz \
   -O /stereoseq/wes/vcf/%.MT.filtered.vcf; \


bcftools view -f PASS /stereoseq/wes/vcf/%.MT.filtered.vcf -Oz -o /stereoseq/wes/vcf/%.MT.filtered.PASS.vcf; \


aws s3 cp /stereoseq/wes/vcf/%.MT.filtered.vcf s3://crm.tumorstudy.mamduh/WES/somatic_variants/GATK4_Mutect2_MT_variants/%.MT.filtered.vcf; \
aws s3 cp /stereoseq/wes/vcf/%.MT.filtered.vcf.idx s3://crm.tumorstudy.mamduh/WES/somatic_variants/GATK4_Mutect2_MT_variants/%.MT.filtered.vcf.idx; \
aws s3 cp /stereoseq/wes/vcf/%.MT.filtered.PASS.vcf s3://crm.tumorstudy.mamduh/WES/somatic_variants/GATK4_Mutect2_MT_variants/%.MT.filtered.PASS.vcf; \
aws s3 cp /stereoseq/wes/vcf/%.MT.filtered.PASS.vcf.idx s3://crm.tumorstudy.mamduh/WES/somatic_variants/GATK4_Mutect2_MT_variants/%.MT.filtered.PASS.vcf.idx; \

rm /stereoseq/wes/bam/*;
rm /stereoseq/wes/vcf/*;

";
