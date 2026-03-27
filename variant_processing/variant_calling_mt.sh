## MITOCHONDRIAL VARIANT CALLING W/ MUTECT2

cd /stereoseq/wes/SD0403/;

#### hs37d5
## filter for MT reads
samtools view -b /stereoseq/wes/SD0403/bam/SD0403_TUM.recalibrated.bam MT > /stereoseq/wes/SD0403/bam/SD0403_TUM.recalibrated.MT.bam
samtools index /stereoseq/wes/SD0403/bam/SD0403_TUM.recalibrated.MT.bam

## run mutect2 on mitochondrial mode
~/tools/gatk-4.6.2.0/gatk Mutect2 \
   -R /stereoseq/reference/hs37d5.fa.gz \
   -I /stereoseq/wes/SD0403/bam/SD0403_TUM.recalibrated.MT.bam \
   --mitochondria-mode true \
   -O /stereoseq/wes/SD0403/vcf/SD0403.wes_tum.MT.vcf

## filtering
~/tools/gatk-4.6.2.0/gatk FilterMutectCalls \
   -V /stereoseq/wes/SD0403/vcf/SD0403.wes_tum.MT.vcf \
   -R /stereoseq/reference/hs37d5.fa.gz \
   -O /stereoseq/wes/SD0403/vcf/SD0403.wes_tum.MT.filtered.vcf

### filter for PASS
bcftools view -f PASS /stereoseq/wes/SD0403/vcf/SD0403.wes_tum.MT.filtered.vcf -Oz -o /stereoseq/wes/SD0403/vcf/SD0403.wes_tum.MT.filtered1.vcf 

### grch38
## filter for MT reads
samtools view -b /stereoseq/wes/SD0403/bam/SD0403.sorted.RG-added.merged.dedup.recal.bam MT > /stereoseq/wes/SD0403/bam/SD0403_TUM.recalibrated.hg38.MT.bam
samtools index /stereoseq/wes/SD0403/bam/SD0403_TUM.recalibrated.hg38.MT.bam

## run mutect2 on mitochondrial mode
~/tools/gatk-4.6.2.0/gatk Mutect2 \
   -R /stereoseq/reference/GRCh38.109.fa \
   -I /stereoseq/wes/SD0403/bam/SD0403_TUM.recalibrated.hg38.MT.bam \
   --mitochondria-mode true \
   -O /stereoseq/wes/SD0403/vcf/SD0403.wes_tum.hg38.MT.vcf

## filtering
~/tools/gatk-4.6.2.0/gatk FilterMutectCalls \
   -V /stereoseq/wes/SD0403/vcf/SD0403.wes_tum.hg38.MT.vcf \
   -R /stereoseq/reference/GRCh38.109.fa \
   -O /stereoseq/wes/SD0403/vcf/SD0403.wes_tum.hg38.MT.filtered.vcf

### filter for PASS
bcftools view -f PASS /stereoseq/wes/SD0403/vcf/SD0403.wes_tum.hg38.MT.filtered.vcf -Oz -o /stereoseq/wes/SD0403/vcf/SD0403.wes_tum.hg38.MT.filtered1.vcf 

######################################################################################
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
