TUM_ID=$1
muttype=$2; ## should be fusion
path=$(echo /stereoseq/all_samples/neoantigens/fusion);

touch $path/logs/$TUM_ID.$muttype.log.txt; ## create log file
linenumber=$(wc -l $path/logs/$TUM_ID.$muttype.log.txt | awk '{print $1}'); ## for resuming job; depracated

## write number of transcripts in each bam into file
echo neopeptide hla chr.pos epitope count > $path/$TUM_ID/counts/$TUM_ID.$muttype.counts.tsv;
while read line; do
    ## define new variant
    neopeptide=$(echo $line | awk '{print $4}');
    hla=$(echo $line | awk '{print $6}');
    chr=$(echo $line | awk '{print $16}');
    start=$(echo $line | awk '{print $17}');
    stop=$(echo $line | awk '{print $18}');
    fullpos=$(echo $chr\:$start\-$stop); ## extract full position
    epitope=$(echo $line | awk '{print $8}');
    echo $neopeptide $hla chr$fullpos $epitope >> $path/logs/$TUM_ID.$muttype.log.txt;

    ## extract transcripts 
    samtools view -b /stereoseq/all_samples/bams/$TUM_ID/$TUM_ID.recal.RGA.bam $fullpos > \
    $path/$TUM_ID/bams/$TUM_ID.$neopeptide.$hla.chr$fullpos.$epitope.$muttype.bam
    
    ## count number of transcripts in each corresponding bam and write into file
    count=$(samtools view $path/$TUM_ID/bams/$TUM_ID.$neopeptide.$hla.chr$fullpos.$epitope.$muttype.bam | wc -l | awk '{print $1}');
    echo $neopeptide $hla chr$fullpos $epitope $count >> $path/$TUM_ID/counts/$TUM_ID.$muttype.counts.tsv; # writing counts to file
done < <(grep $TUM_ID $path/stereo_fusion_neopeptides_sep.tsv) ## input defined by $muttype, default should be 'all' to extract all somatic.filtered1.funcotated.vcf variants

touch $path/logs/$TUM_ID.$muttype.extract_reads.success;
