TUM_ID=$1;
muttype=$2;
path=$(echo /stereoseq/all_samples/neoantigens/fusion);

# for file in $(ls $path/bams/$muttype/*.$muttype.bam); do
#     touch $path/coords/$TUM_ID.neoantigens.$muttype.coords.txt;
#     samtools view $file | awk '{print $1,$(NF-1),$NF}' >> $path/coords/$TUM_ID.neoantigens.$muttype.coords.txt;
# done

touch $path/$TUM_ID/coords/$TUM_ID.$muttype.neoantigens.confirmed.coords.txt;
while read line; do 
    ## extract information from FN
	neopeptide=$(echo $line | awk '{print $4}');
    gene=$(echo $line | awk '{print $5}');
    sequence=$(echo $line | awk '{print $6}');
    hla=$(echo $line | awk '{print $7}');
    epitope=$(echo $line | awk '{print $9}');
    variantType=$(echo $line | awk '{print $10}');
    chr=$(echo $line | awk '{print $17}');
    start=$(echo $line | awk '{print $18}');
    stop=$(echo $line | awk '{print $19}');
    fullpos=$(echo $chr\:$start\-$stop);
    transcript=$(echo $line | awk '{print $20}');
    
    ## define corresponding bam file
    file=$path/$TUM_ID/bams/$TUM_ID.$neopeptide.$hla.chr$fullpos.$epitope.$muttype.bam

    ## get sliding window of 20bp each
    echo "$sequence" | awk -v k=20 '{
    for (i=1; i<=length($0)-k+1; i++)
        print substr($0,i,k)
    }' > $path/$TUM_ID/coords/$gene.windows.txt

    samtools view "$file" | grep -F -f $path/$TUM_ID/coords/$gene.windows.txt | while read line2; do

        transcriptID=$(echo "$line2" | awk '{print $1}')
        readseq=$(echo "$line2" | awk '{print $10}')
        x=$(echo "$line2" | awk '{print $(NF-1)}')
        y=$(echo "$line2" | awk '{print $NF}')

        echo "$neopeptide,$gene,$sequence,$hla,$epitope,$variantType,$chr,$start,$stop,$transcript,\
              $transcriptID,$readseq,$x,$y" >> "$path/$TUM_ID/coords/$TUM_ID.$muttype.neoantigens.confirmed.coords.txt"
    done #< <(samtools view $file | grep $subseq)
done < <(grep $TUM_ID $path/fusion_neoantigens_confirmed_sequence.tsv)
