TUM_ID=$1
muttype=$2;
in_file=$3
path=$(echo /stereoseq/all_samples/mutations/$TUM_ID);

# for file in $(ls $path/bams/$muttype/*.$muttype.bam); do
#     touch $path/coords/$TUM_ID.somatic.$muttype.coords.txt;
#     samtools view $file | awk '{print $1,$(NF-1),$NF}' >> $path/coords/$TUM_ID.somatic.$muttype.coords.txt;
# done

> $path/coords/$TUM_ID.somatic.$muttype.coords.txt; ## create coords file
while read line; do 
    ## extract information from VCF
	chr=$(echo $line | awk '{print $1}');
    pos=$(echo $line | awk '{print $2}');
    ref=$(echo $line | awk '{print $4}');
    alt=$(echo $line | awk '{print $5}');
    gene=$(echo $line | awk '{print $8}' | awk -F'[][]' '{match($2, /^[^|]+/, a); print a[0]}');
    variantClassification=$(echo $line | awk '{print $8}' | cut -d '|' -f 7);
    variantType=$(echo $line | awk '{print $8}' | cut -d '|' -f 9);
    annotationTranscript=$(echo $line | awk '{print $8}' | cut -d '|' -f 14);
    transcriptStrand=$(echo $line | awk '{print $8}' | cut -d '|' -f 15);
    genomeChange=$(echo $line | awk '{print $8}' | cut -d '|' -f 13);
    cDnaChange=$(echo $line | awk '{print $8}' | cut -d '|' -f 18);
    codonChange=$(echo $line | awk '{print $8}' | cut -d '|' -f 19);
    proteinChange=$(echo $line | awk '{print $8}' | cut -d '|' -f 20);
    context=$(echo $line | cut -d '|' -f 22);

    ## define corresponding bam file
    file=$path/bams/all/$TUM_ID.somatic.chr$chr.$pos.all.bam;
    
    while read line2; do
        ## extract information from corresponding bam file
        transcriptID=$(echo $line2 | awk '{print $1}')
        x=$(echo $line2 | awk '{print $(NF-1)}') ## X coordinate
        y=$(echo $line2 | awk '{print $NF}') ## Y coordinate

        ## write VCF + bam information into file
        echo $chr,$pos,$ref,$alt,$gene,$context,$variantClassification,\
        $variantType,$annotationTranscript,$transcriptStrand,$genomeChange,$cDnaChange,$codonChange,$proteinChange,\
        $transcriptID,$x,$y >> $path/coords/$TUM_ID.somatic.$muttype.coords.txt;
    done < <(samtools view $file)
done < /stereoseq/all_samples/vcf/$TUM_ID/$TUM_ID.$in_file.vcf 
## input file defined as var to include only subset of variants, e.g. filtered7, MT, etc.
## input VCF should have NO headers