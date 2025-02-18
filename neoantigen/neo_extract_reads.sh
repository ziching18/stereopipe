TUM_ID=$1;
path=$(echo /stereoseq/all_samples/neoantigens/${TUM_ID});
touch $path/logs/$TUM_ID.neoantigens.log.txt;
linenumber=$(wc -l $path/logs/$TUM_ID.neoantigens.log.txt | awk '{print $1}');

echo chr pos ref alt gene chr.pos.ref.alt.gene count > $path/counts/$TUM_ID.neoantigens.vcf.transcript.gene.counts.tsv;
while read line; do
    rm $path/subset.bam;
    rm $path/skrip.js;
    chr=$(echo $line | awk '{print $1}');
    pos=$(echo $line | awk '{print $2}');
    fullpos=$(echo $chr\:$(($pos-200))\-$(($pos+200)));
    echo $fullpos >> $path/${TUM_ID}.neoantigens.log.txt;
    ref=$(echo $line | awk '{print $5}');
    alt=$(echo $line | awk '{print $6}' | cut -c1-1);
    gene=$(echo $line | awk '{print $7}');
    echo "final String contig= \"$chr\";" >> $path/skrip.js;
    echo "final int mutpos = $pos;" >> $path/skrip.js;
    echo "final char mutbase='$alt';" >> $path/skrip.js;
    cat $path/scriptbody.js >> $path/skrip.js;
    samtools view -b /stereoseq/all_samples/bams/${TUM_ID}/${TUM_ID}.recal.RGA.bam $fullpos > $path/subset.bam;
    java -jar ~/tools/jvarkit/dist/jvarkit.jar samjdk \
    -f $path/skrip.js $path/subset.bam \
    -o $path/bams/${TUM_ID}.neoantigens.chr$chr.$pos.bam --samoutputformat BAM;
    count=$(samtools view $path/bams/$muttype/$TUM_ID.somatic.chr$chr.$pos.$muttype.bam | wc -l | awk '{print $1}');
    echo $chr $pos $ref $alt $gene chr$chr.$pos.$ref.$alt.$gene $count >> $path/counts/$TUM_ID.neoantigens.vcf.transcript.gene.counts.tsv; # writing counts to file
done < $path/${TUM_ID}.neoantigens.filtered.full.vcf.tsv
#<(tail -n "+$linenumber" $path/${TUM_ID}.neoantigens.${CLASS}.filtered.full.vcf.tsv)
