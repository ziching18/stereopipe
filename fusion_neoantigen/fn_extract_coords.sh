TUM_ID=$1;
muttype=$2;
path=$(echo /stereoseq/all_samples/neoantigens/fusion);

# for file in $(ls $path/bams/$muttype/*.$muttype.bam); do
#     touch $path/coords/$TUM_ID.neoantigens.$muttype.coords.txt;
#     samtools view $file | awk '{print $1,$(NF-1),$NF}' >> $path/coords/$TUM_ID.neoantigens.$muttype.coords.txt;
# done

touch $path/$TUM_ID/coords/$TUM_ID.$muttype.neoantigens.coords.txt;
while read line; do 
    ## extract information from FN
	neopeptide=$(echo $line | awk '{print $4}');
    hla=$(echo $line | awk '{print $6}');
    gene=$(echo $line | awk '{print $5}');
    epitope=$(echo $line | awk '{print $8}');
    variantType=$(echo $line | awk '{print $9}');
    chr=$(echo $line | awk '{print $16}');
    start=$(echo $line | awk '{print $17}');
    stop=$(echo $line | awk '{print $18}');
    fullpos=$(echo $chr\:$start\-$stop);
    transcript=$(echo $line | awk '{print $19}');
    
    ## define corresponding bam file
    file=$path/$TUM_ID/bams/$TUM_ID.$neopeptide.$hla.chr$fullpos.$epitope.$muttype.bam
    
    while read line2; do
        ## extract information from corresponding bam file
        transcriptID=$(echo $line2 | awk '{print $1}')
        x=$(echo $line2 | awk '{print $(NF-1)}')
        y=$(echo $line2 | awk '{print $NF}')

        ## write FN + bam information into file
        echo $neopeptide,$hla,$gene,$epitope,$variantType,$chr,$start,$stop,$transcript,\
        $transcriptID,$x,$y >> $path/$TUM_ID/coords/$TUM_ID.$muttype.neoantigens.coords.txt;
    done < <(samtools view $file)
done < <(grep $TUM_ID $path/stereo_fusion_neopeptides_sep.tsv)
