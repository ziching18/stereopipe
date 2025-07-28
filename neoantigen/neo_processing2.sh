TUM_ID=$1;
# rows=$(grep "#" /stereoseq/all_samples/vcf/${TUM_ID}/${TUM_ID}.somatic.funcotated.vcf | wc -l);

# python3 /stereoseq/code/stereopipe/variant_calling/vcf_filter_variants.py -s ${TUM_ID} \
# --skiprows $rows --in_dir "/stereoseq/all_samples/vcf/${TUM_ID}" \
# --out_dir "/stereoseq/all_samples/mutations/${TUM_ID}/variants";

# bash /stereoseq/code/stereopipe/neoantigen/neo_extract_vcf.sh ${TUM_ID};

bash /stereoseq/code/stereopipe/neoantigen/neo_extract_coords.sh ${TUM_ID} filtered;

python3 /stereoseq/code/stereopipe/neoantigen/neo_raw_round.py -s ${TUM_ID} \
--bin_size 100 --dir "/stereoseq/all_samples/neoantigens/${TUM_ID}" --muttype filtered;