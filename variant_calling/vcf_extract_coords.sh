TUM_ID=$1
muttype=$2;
path=$(echo /stereoseq/all_samples/mutations/$TUM_ID);

# for file in $(ls $path/bams/$muttype/*.$muttype.bam); do
#     touch $path/coords/$TUM_ID.somatic.$muttype.coords.txt;
#     samtools view $file | awk '{print $1,$(NF-1),$NF}' >> $path/coords/$TUM_ID.somatic.$muttype.coords.txt;
# done

> $path/coords/$TUM_ID.somatic.$muttype.coords.txt;
while read line; do 
	chr=$(echo $line | awk '{print $1}');
    pos=$(echo $line | awk '{print $2}');
    ref=$(echo $line | awk '{print $4}');
    alt=$(echo $line | awk '{print $5}');
    gene=$(echo $line | awk '{print $8}' | awk -F'[][]' '{match($2, /^[^|]+/, a); print a[0]}');
    variantClassification=$(echo $line | awk '{print $8}' | cut -d '|' -f 7);
    variantType=$(echo $line | awk '{print $8}' | cut -d '|' -f 9);
    genomeChange=$(echo $line | awk '{print $8}' | cut -d '|' -f 13);
    cDnaChange=$(echo $line | awk '{print $8}' | cut -d '|' -f 18);
    codonChange=$(echo $line | awk '{print $8}' | cut -d '|' -f 19);
    proteinChange=$(echo $line | awk '{print $8}' | cut -d '|' -f 20);
    context=$(echo $line | cut -d '|' -f 22);
    file=$path/bams/all/$TUM_ID.somatic.chr$chr.$pos.all.bam;
    
    while read line2; do
        transcriptID=$(echo $line2 | awk '{print $1}')
        x=$(echo $line2 | awk '{print $(NF-1)}')
        y=$(echo $line2 | awk '{print $NF}')
        echo $chr,$pos,$ref,$alt,$gene,$context,$variantClassification,\
        $variantType,$genomeChange,$cDnaChange,$codonChange,$proteinChange,\
        $transcriptID,$x,$y >> $path/coords/$TUM_ID.somatic.$muttype.coords.txt;
    done < <(samtools view $file)
done < /stereoseq/all_samples/vcf/$TUM_ID/$TUM_ID.somatic.filtered5.funcotated.vcf