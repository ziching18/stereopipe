TUM_ID=$1;

bash /stereoseq/code/stereopipe/neoantigen/neo_extract_vcf.sh ${TUM_ID};

python3 /stereoseq/code/stereopipe/neoantigen/neo_filter_variants.py -s ${TUM_ID} \
--in_dir "/stereoseq/all_samples/neoantigens/${TUM_ID}" \
--out_dir "/stereoseq/all_samples/neoantigens/${TUM_ID}/variants";

for muttype in all CT GA AT TCA YTCA RTCA; do
    if [ ! -f /stereoseq/all_samples/neoantigens/$TUM_ID/logs/$TUM_ID.$muttype.extract_reads.success ]; then
        bash /stereoseq/code/stereopipe/neoantigen/neo_extract_reads.sh ${TUM_ID} ${muttype};
    else
        echo "Reads already extracted sis";
    fi;
done;

for muttype in all CT GA AT TCA YTCA RTCA; do
    if [ ! -f /stereoseq/all_samples/neoantigens/$TUM_ID/coords/$TUM_ID.somatic.$muttype.coords.txt ]; then
        bash /stereoseq/code/stereopipe/neoantigen/neo_extract_coords.sh ${TUM_ID} ${muttype};
    else
        echo "Coords already extracted sis";
    fi;
done;

python3 /stereoseq/code/stereopipe/neoantigen/neo_round_counts.py -s ${TUM_ID} \
--bin_size 100 --dir "/stereoseq/all_samples/neoantigens/${TUM_ID}";

python3 /stereoseq/code/stereopipe/neoantigen/neo_normlog.py -s ${TUM_ID} \
--bin_size 100 --dir "/stereoseq/all_samples/neoantigens/${TUM_ID}";