TUM_ID=$1;

lines1=$(wc -l < /stereoseq/all_samples/neoantigens/$TUM_ID/$TUM_ID_stereo_TUM.filtered.tsv)
lines2=$(($lines1 - 1))
echo $(head -1 /stereoseq/all_samples/neoantigens/$TUM_ID/${TUM_ID}_stereo_TUM.filtered.tsv) \
$(grep "#" /stereoseq/all_samples/vcf/$TUM_ID/$TUM_ID.somatic.filtered7.funcotated.vcf | tail -1) \
> /stereoseq/all_samples/neoantigens/$TUM_ID/$TUM_ID.neoantigens.all.vcf.tsv;

while read line; do
	chr=$(echo $line | awk '{print $1}');
	var_stop=$(echo $line | awk '{print $3}');
	transcript=$(echo $line | awk '{print $6}');
	peptide=$(echo $line | awk '{print $21}');
	vcfLine=$(awk -v stop="$var_stop" '$2 ~ stop {print}' /stereoseq/all_samples/vcf/$TUM_ID/$TUM_ID.somatic.filtered7.funcotated.vcf);
    echo $vcfLine $line >> /stereoseq/all_samples/neoantigens/$TUM_ID/$TUM_ID.neoantigens.all.vcf.tsv;
done < <(tail -$lines2 /stereoseq/all_samples/neoantigens/$TUM_ID/${TUM_ID}_stereo_TUM.filtered.tsv)
