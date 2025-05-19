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
    ref=$(echo $line | awk '{print $4}');
    alt=$(echo $line | awk '{print $5}');
    gene=$(echo $line | awk '{print $8}' | awk -F'[][]' '{match($2, /^[^|]+/, a); print a[0]}')
    file=$path/bams/all/$TUM_ID.somatic.chr$chr.$pos.all.bam
    while read line2; do
        echo $chr $pos $ref $alt $gene $(echo $line2 | awk '{print $1,$(NF-1),$NF}') >> $path/coords/$TUM_ID.somatic.$muttype.coords.txt;
    done < <(samtools view $file)
done < /stereoseq/all_samples/vcf/$TUM_ID/$TUM_ID.somatic.filtered5.funcotated.vcf