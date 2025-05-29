# NEOANTIGEN PREDICTION
# while read line; do 
#	slx=$(echo $line | awk '{print $1}'); 
#	id=$(echo $line | awk '{print $2}'); 
#	hla=$(grep $id /stereoseq/neoantigens/data/neoantigens_hla_alleles.txt | awk '{print}');
#	sudo docker exec neo bash -c bash /stereoseq/code/pvacseq_neoantigens_prediction.sh $slx $id $hla; 
#done < /stereoseq/code/samples.txt

TUM_ID=$1;
HLA=$2; 

echo "running..." ${TUM_ID} ${HLA};
pvacseq run \
/vcf/${TUM_ID}.somatic.funcotated.vep.vcf.gz \
${TUM_ID}_stereo_TUM \
${HLA} \
MHCflurry MHCflurryEL NNalign NetMHC NetMHCIIpan NetMHCIIpanEL NetMHCpan NetMHCpanEL \
/output \
-e1 8,9,10,11 \
-e2 15 \
--normal-sample-name ${TUM_ID}_MN \
-p /vcf/${TUM_ID}.phased.vep.vcf.gz \
--iedb-install-directory /opt/iedb ;
echo "finished ${TUM_ID}"
