TUM_ID=$1
muttype=$2;
path=$(echo /stereoseq/all_samples/mutations/$TUM_ID);
touch $path/logs/$TUM_ID.$muttype.log.txt;
linenumber=$(wc -l $path/logs/$muttype.log.txt | awk '{print $1}');

while read line; do
    rm $path/$muttype.subset.bam;
    rm $path/$muttype.skrip.js;
    chr=$(echo $line | awk '{print $2}');
    pos=$(echo $line | awk '{print $3}');
    fullpos=$(echo $chr\:$(($pos-500))\-$(($pos+500)));
    echo $chr $fullpos >> $path/logs/$muttype.log.txt;
    alt=$(echo $line | awk '{print $6}');
    echo "final String contig= \"$chr\";" >> $path/bams/$muttype.skrip.js;
    echo "final int mutpos = $pos;" >> $path/bams/$muttype.skrip.js;
    echo "final char mutbase='$alt';" >> $path/bams/$muttype.skrip.js;
    cat $path/scriptbody.js >> $path/bams/$muttype.skrip.js;
    samtools view -b /stereoseq/all_samples/bams/$TUM_ID/$TUM_ID.recal.RGA.bam $fullpos > $path/bams/$muttype.subset.bam;
    java -jar ~/tools/jvarkit/dist/jvarkit.jar samjdk \
    -f $path/bams/$muttype.skrip.js $path/bams/$muttype.subset.bam \
    -o $path/bams/$muttype/$TUM_ID.somatic.chr$chr.$pos.$muttype.bam --samoutputformat BAM;
#done < $path/$TUM_ID.somatic.$muttype.tsv
done < <(tail -n "+$linenumber" $path/variants/$TUM_ID.somatic.$muttype.tsv)
