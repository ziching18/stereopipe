TUM_ID=$1
muttype=$2;
path=$(echo /stereoseq/all_samples/mutations/$TUM_ID);

# for file in $(ls $path/bams/$muttype/*.$muttype.bam); do
#     touch $path/coords/$TUM_ID.somatic.$muttype.coords.txt;
#     samtools view $file | awk '{print $1,$(NF-1),$NF}' >> $path/coords/$TUM_ID.somatic.$muttype.coords.txt;
# done

touch $path/coords/$TUM_ID.somatic.$muttype.coords.txt;
while read line; do 
	chr=$(echo $line | awk '{print $1}');
    pos=$(echo $line | awk '{print $2}');
    file=$($path/bams/all/$TUM_ID.somatic.chr$chr.$pos.all.bam)
    samtools view $file | awk '{print $1,$(NF-1),$NF}' >> $path/coords/$TUM_ID.somatic.$muttype.coords.txt;
done < $TUM_ID.somatic.filtered5.funcotated.vcf