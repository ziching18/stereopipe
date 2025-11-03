# NEOANTIGEN PREDICTION
## resume pvacseq neoantigen prediction with created docker instance
## to run in docker environment

TUM_ID=$1;

hla=$(grep ${TUM_ID} data/all_neoantigens_hla_alleles.txt | awk '{print $NF}');

sudo docker exec ${TUM_ID}_neo bash -c "bash /code/neo_pvacseq_prediction.sh ${TUM_ID} $hla";