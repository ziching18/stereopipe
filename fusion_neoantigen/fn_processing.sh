TUM_ID=$1;

bash stereopipe/fusion_neoantigen/fn_extract_reads.sh ${TUM_ID} fusion;

bash stereopipe/fusion_neoantigen/fn_extract_coords.sh ${TUM_ID} fusion;

python3 stereopipe/fusion_neoantigen/fn_raw_round.py -s ${TUM_ID} \
--bin_size 100 --dir "/stereoseq/all_samples/neoantigens/fusion/${TUM_ID}" --muttype fusion;