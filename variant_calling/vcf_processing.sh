TUM_ID=$1
rows=$(grep "#" /stereoseq/all_samples/vcf/${TUM_ID}/${TUM_ID}.somatic.funcotated.vcf | wc -l)

python3 stereopipe/variant_calling/vcf_filter_variants.py -s ${TUM_ID} \
--skiprows $rows --in_dir "/stereoseq/all_samples/vcf/${TUM_ID}/" \
--out_dir "/stereoseq/all_samples/mutations/${TUM_ID}/variants/";

for muttype in all CT GA AT TCA YTCA RTCA TP53 PIK3CA GATA3 MAP3K1 KMT2C PTEN CBFB CDH1 AKT1 NF1; do
    bash stereopipe/variant_calling/vcf_extract_reads.sh ${TUM_ID} ${muttype};
done;