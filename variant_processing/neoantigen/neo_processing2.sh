TUM_ID=$1;
rows=$(grep "#" /stereoseq/all_samples/vcf/${TUM_ID}/${TUM_ID}.somatic.funcotated.vcf | wc -l);

## 1. Convert TSV to VCF
bash stereopipe/variant_processing/neoantigen/neo_extract_vcf.sh ${TUM_ID};

## No need to re-extract bams (performed in VCF processing pipeline)

bash stereopipe/variant_processing/neoantigen/neo_extract_coords.sh ${TUM_ID} filtered;

python3 stereopipe/variant_processing/neoantigen/neo_raw_round.py -s ${TUM_ID} \
--bin_size 100 --dir "/stereoseq/all_samples/neoantigens/${TUM_ID}" --muttype filtered;