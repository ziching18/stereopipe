TUM_ID=$1;
rows=$(grep "#" /stereoseq/all_samples/vcf/${TUM_ID}/${TUM_ID}.somatic.funcotated.vcf | wc -l);

bash /stereoseq/code/stereopipe/neoantigens/neo_extract_vcf.sh ${TUM_ID};

python3 /stereoseq/code/stereopipe/variant_calling/neo_filter_variants.py -s ${TUM_ID} \
--in_dir "/stereoseq/all_samples/neoantigens/${TUM_ID}" \
--out_dir "/stereoseq/all_samples/neoantigens/${TUM_ID}/variants";

for muttype in all CT GA AT TCA YTCA RTCA TP53 PIK3CA GATA3 MAP3K1 KMT2C PTEN CBFB CDH1 AKT1 NF1; do
    if [ ! -f /stereoseq/all_samples/mutations/$TUM_ID/logs/$TUM_ID.$muttype.extract_reads.success ]; then
        bash /stereoseq/code/stereopipe/variant_calling/neo_extract_reads.sh ${TUM_ID} ${muttype};
    else
        echo "Reads already extracted sis";
    fi;
done;

for muttype in all CT GA AT TCA YTCA RTCA TP53 PIK3CA GATA3 MAP3K1 KMT2C PTEN CBFB CDH1 AKT1 NF1; do
    if [ ! -f /stereoseq/all_samples/mutations/$TUM_ID/coords/$TUM_ID.somatic.$muttype.coords.txt ]; then
        bash /stereoseq/code/stereopipe/variant_calling/neo_extract_coords.sh ${TUM_ID} ${muttype};
    else
        echo "Coords already extracted sis";
    fi;
done;

python3 /stereoseq/code/stereopipe/variant_calling/neo_round_counts.py -s ${TUM_ID} \
--bin_size 100 --dir "/stereoseq/all_samples/mutations/${TUM_ID}";

python3 /stereoseq/code/stereopipe/variant_calling/neo_normlog.py -s ${TUM_ID} \
--bin_size 100 --dir "/stereoseq/all_samples/mutations/${TUM_ID}";