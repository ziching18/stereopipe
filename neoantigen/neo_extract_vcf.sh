TUM_ID=$1;

cd /stereoseq/all_samples/neoantigens/$TUM_ID;

## remove spaces in the header
sed -i 's/ /_/g' /stereoseq/all_samples/neoantigens/$TUM_ID/${TUM_ID}_stereo_TUM.filtered.tsv;
# | head -1 \
> /stereoseq/all_samples/neoantigens/$TUM_ID/$TUM_ID.neoantigens.all.vcf.tsv;

lines1=$(wc -l < /stereoseq/all_samples/neoantigens/$TUM_ID/${TUM_ID}_stereo_TUM.filtered.tsv)
lines2=$(($lines1 - 1))
while read line; do
	chr=$(echo $line | awk '{print $1}');
	var_stop=$(echo $line | awk '{print $3}');
	transcript=$(echo $line | awk '{print $6}');
	peptide=$(echo $line | awk '{print $21}');
	awk -v stop="$var_stop" '$2 ~ stop {print}' /stereoseq/all_samples/vcf/$TUM_ID/$TUM_ID.somatic.filtered7.funcotated.vcf \
    >> /stereoseq/all_samples/neoantigens/$TUM_ID/$TUM_ID.neoantigens.all.vcf.tsv;
done < <(tail -$lines2 /stereoseq/all_samples/neoantigens/$TUM_ID/${TUM_ID}_stereo_TUM.filtered.tsv)
