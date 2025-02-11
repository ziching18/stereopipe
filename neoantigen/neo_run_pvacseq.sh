TUM_ID=$1;

mkdir /stereoseq/neoantigens/output_${TUM_ID};

hla=$(grep ${TUM_ID} /stereoseq/neoantigens/data/neoantigens_hla_alleles.txt | awk '{print $NF}');
sudo docker rm --force ${TUM_ID}_neo;
sudo docker run \
--name ${TUM_ID}_neo \
-v /stereoseq/all_samples/vcf:/vcf \
-v /stereoseq/neoantigens/data:/pvactools_stereoseq_data \
-v /stereoseq/neoantigens/output_${TUM_ID}:/pvactools_stereoseq_output \
-v /home/ubuntu/tools/VEP_plugins:/VEP_plugins \
-v /stereoseq/code:/code \
-dit griffithlab/pvactools ;
sudo docker exec ${TUM_ID}_neo bash -c "bash /code/pvacseq_neoantigens_prediction.sh ${TUM_ID} $hla"; 
