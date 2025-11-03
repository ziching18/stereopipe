TUM_ID=$1;
FORCE=$2;

## read FilterMutectCalls + Funcotator VCF
rows=$(grep "#" /stereoseq/all_samples/vcf/${TUM_ID}/${TUM_ID}.somatic.filtered1.funcotated.vcf | wc -l);

## 1. Extract VCF information (all): vcf_extract_info.py
python3 stereopipe/variant_processing/somatic_mutation/vcf_extract_info.py -s ${TUM_ID} \
--skiprows $rows --in_dir "/stereoseq/all_samples/vcf/${TUM_ID}" \
--out_dir "/stereoseq/all_samples/mutations/${TUM_ID}/variants" $FORCE;

## 2. Extract reads from corresponding bams (all): vcf_extract_reads.sh
bash stereopipe/variant_processing/somatic_mutation/vcf_extract_reads.sh ${TUM_ID} all;

## 3. Extract coordinates from corresponding bams (specific subsets), input VCF should have NO header: vcf_extract_coords.sh
bash stereopipe/variant_processing/somatic_mutation/vcf_extract_coords.sh ${TUM_ID} filtered somatic.filtered5;

## 4. Round coordinates to nearest bin_size: vcf_raw_round.py
python3 stereopipe/variant_processing/somatic_mutation/vcf_raw_round.py -s ${TUM_ID} \
--bin_size 100 --dir "/stereoseq/all_samples/mutations/${TUM_ID}" --muttype filtered ;