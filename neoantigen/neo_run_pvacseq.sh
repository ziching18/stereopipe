TUM_ID=$1;

#mkdir /mnt/c/Users/staffgenomics/wsl_workdir/ziching/neo/code/out/output_${TUM_ID};

hla=$(grep ${TUM_ID} data/all_neoantigens_hla_alleles.txt | awk '{print $NF}');
echo ${TUM_ID} $hla;
#sudo docker rm --force ${TUM_ID}_neo;
sudo docker run \
--name ${TUM_ID}_neo \
-v /home/ziching/neo/vcf:/vcf \
-v /home/ziching/neo/data:/data \
-v /home/ziching/neo/out/output_${TUM_ID}:/output \
-v /home/ziching/neo/VEP_plugins:/VEP_plugins \
-v /home/ziching/neo/code/stereopipe/neoantigen:/code \
-dit griffithlab/pvactools \
--restart=unless-stopped ;
sudo docker exec ${TUM_ID}_neo bash -c "bash /code/neo_pvacseq_prediction.sh ${TUM_ID} $hla"; 
