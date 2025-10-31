TUM_ID=$1
muttype=$2;
path=$(echo /stereoseq/all_samples/mutations/$TUM_ID);

touch $path/logs/$TUM_ID.$muttype.log.txt; ## create log file
linenumber=$(wc -l $path/logs/$TUM_ID.$muttype.log.txt | awk '{print $1}'); ## for resuming job; depracated

## write number of transcripts in each bam into file
echo chr pos ref alt gene chr.pos.ref.alt.gene count > $path/counts/$TUM_ID.somatic.vcf.transcript.$muttype.gene.counts.tsv;
while read line; do
    ## remove previous temp bam & script files
    rm $path/bams/$muttype.subset.bam;
    rm $path/bams/$muttype.skrip.js;

    ## define new variant
    chr=$(echo $line | awk '{print $2}');
    pos=$(echo $line | awk '{print $3}');
    fullpos=$(echo $chr\:$(($pos-500))\-$(($pos+500))); ## extract 500 positions before and after base
    echo $chr $fullpos >> $path/logs/$TUM_ID.$muttype.log.txt;
    ref=$(echo $line | awk '{print $5}');
    alt=$(echo $line | awk '{print $6}');
    gene=$(echo $line | awk '{print $7}');

    ## define variables for samjdk into temporary script file: 'skrip.js'
    echo "final String contig= \"$chr\";" >> $path/bams/$muttype.skrip.js;
    echo "final int mutpos = $pos;" >> $path/bams/$muttype.skrip.js;
    echo "final char mutbase='$alt';" >> $path/bams/$muttype.skrip.js;
    cat /stereoseq/all_samples/mutations/scriptbody.js >> $path/bams/$muttype.skrip.js;

    ## extract transcripts into temporary bam file: 'subset.bam'
    samtools view -b /stereoseq/all_samples/bams/$TUM_ID/$TUM_ID.recal.RGA.bam $fullpos > $path/bams/$muttype.subset.bam;
    
    ## run samjdk with temporary bam and script files, 
    ## output corresponding bams into bams/all/
    java -jar ~/tools/jvarkit/dist/jvarkit.jar samjdk \
    -f $path/bams/$muttype.skrip.js $path/bams/$muttype.subset.bam \
    -o $path/bams/$muttype/$TUM_ID.somatic.chr$chr.$pos.$muttype.bam --samoutputformat BAM;
    
    ## count number of transcripts in each corresponding bam and write into file
    count=$(samtools view $path/bams/$muttype/$TUM_ID.somatic.chr$chr.$pos.$muttype.bam | wc -l | awk '{print $1}');
    echo $chr $pos $ref $alt $gene chr$chr.$pos.$ref.$alt.$gene $count >> $path/counts/$TUM_ID.somatic.vcf.transcript.$muttype.gene.counts.tsv; # writing counts to file
#done < <(tail -n "+$linenumber" $path/variants/$TUM_ID.somatic.$muttype.tsv)
done < $path/variants/$TUM_ID.somatic.$muttype.tsv ## input defined by $muttype, default should be 'all' to extract all somatic.filtered1.funcotated.vcf variants
touch $path/logs/$TUM_ID.$muttype.extract_reads.success;
