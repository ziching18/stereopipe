TUM_ID=$1;

path=$(echo /stereoseq/all_samples/neoantigens/$TUM_ID);

## remove spaces in the header
sed -i 's/ /_/g' $path/${TUM_ID}_stereo_TUM.filtered.tsv;
# | head -1 \

## extract df header from source VCF  and write into file (file creation)
grep "#" /stereoseq/all_samples/vcf/$TUM_ID/$TUM_ID.somatic.filtered7.funcotated.vcf | tail -1 > $path/$TUM_ID.somatic.neoantigens.all.vcf.tsv;

lines1=$(wc -l < $path/${TUM_ID}_stereo_TUM.filtered.tsv)
lines2=$(($lines1 - 1))
while read line; do
	## extract information from neoantigen TSV
	var_stop=$(echo $line | awk '{print $3}');
	transcript=$(echo $line | awk '{print $6}');
	hlaAllele=$(echo $line | awk '{print $17}');
	MTpeptide=$(echo $line | awk '{print $21}');
	WTpeptide=$(echo $line | awk '{print $22}');

	## extract information from corresponding variant, from second to last column
	variant=$(awk -v stop="$var_stop" '$2 ~ stop {print}' /stereoseq/all_samples/vcf/$TUM_ID/$TUM_ID.somatic.filtered7.funcotated.vcf)

	## write all information into file
	echo $variant $transcript $hlaAllele $MTpeptide $WTpeptide >> $path/$TUM_ID.somatic.neoantigens.all.vcf.tsv;
done < <(tail -$lines2 $path/${TUM_ID}_stereo_TUM.filtered.tsv) ## exclude header
