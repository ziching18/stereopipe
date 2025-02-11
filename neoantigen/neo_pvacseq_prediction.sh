## NEOANTIGEN PREDICTION

TUM_ID=$1;
HLA=$2; 

#echo "running..." ${TUM_ID} ${HLA};
pvacseq run \
/vcf/${TUM_ID}/${TUM_ID}.somatic.funcotated.vep.vcf.gz \
${TUM_ID} \
${HLA} \
MHCflurry MHCflurryEL NNalign NetMHC NetMHCIIpan NetMHCIIpanEL NetMHCpan NetMHCpanEL \
/pvactools_stereoseq_output \
-e1 8,9,10,11 \
-e2 15 \
-p /vcf/${TUM_ID}/${TUM_ID}.phased.vep.vcf.gz \
--iedb-install-directory /opt/iedb ;
#echo "finished ${TUM_ID}