TUM_ID=$1

cd /stereoseq/all_samples/vcf/for_sig/;

# Step 0: normalize multiallelics
bcftools norm -m -any /stereoseq/all_samples/vcf/${TUM_ID}/${TUM_ID}.somatic.filtered1.vcf -Oz -o step0.vcf

# Step 1: keep PASS + SNPs only
bcftools view -f PASS -v snps step0.vcf -Oz -o step1.vcf

# Step 2: quality filters
bcftools view -i 'FORMAT/DP>10 && FORMAT/AD[1]>3 && FORMAT/AF[1:0]>=0.1 && FORMAT/AF[1:0]<=0.8 && INFO/TLOD>=6' step1.vcf -Oz -o tmp.vcf

bcftools view -e 'INFO/STR=1' tmp.vcf -Oz -o step2.vcf

# Step 3: remove problematic regions
bcftools view -r 1-22 step2.vcf -Oz -o step3.vcf

# Step 4: optional IG/HLA removal
grep -v -e "\[IGH" -e "\[IGK" -e "\[IGL" -e "\[HLA" step3.vcf > ${TUM_ID}.final.vcf