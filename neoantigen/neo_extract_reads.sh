TUM_ID=$1;
muttype=$2;

path=$(echo /stereoseq/all_samples/neoantigens/$TUM_ID);
touch $path/logs/$TUM_ID.neoantigens.$muttype.log.txt;
linenumber=$(wc -l $path/logs/$TUM_ID.neoantigens.$muttype.log.txt | awk '{print $1}');

echo chr pos ref alt gene chr.pos.ref.alt.gene count > $path/counts/$TUM_ID.neoantigens.vcf.transcript.$muttype.gene.counts.tsv;
while read line; do
    rm $path/bams/$muttype.subset.bam;
    rm $path/bams/$muttype.skrip.js;
    chr=$(echo $line | awk '{print $2}');
    pos=$(echo $line | awk '{print $3}');
    fullpos=$(echo $chr\:$(($pos-500))\-$(($pos+500)));
    echo $chr $fullpos >> $path/logs/$TUM_ID.neoantigens.$muttype.log.txt;
    ref=$(echo $line | awk '{print $5}');
    alt=$(echo $line | awk '{print $6}' | cut -c1-1);
    #gene=$(echo $line | awk '{print $8}' | cut -d '[' -f 2 | cut -d '|' -f 1);
    gene=$(echo $line | awk '{print $7}');
    echo "final String contig= \"$chr\";" >> $path/bams/$muttype.skrip.js;
    echo "final int mutpos = $pos;" >> $path/bams/$muttype.skrip.js;
    echo "final char mutbase='$alt';" >> $path/bams/$muttype.skrip.js;
    cat /stereoseq/all_samples/neoantigens/scriptbody.js >> $path/bams/$muttype.skrip.js;
    samtools view -b /stereoseq/all_samples/bams/$TUM_ID/$TUM_ID.recal.RGA.bam $fullpos > $path/bams/$muttype.subset.bam;
    java -jar ~/tools/jvarkit/dist/jvarkit.jar samjdk \
    -f $path/bams/$muttype.skrip.js $path/bams/$muttype.subset.bam \
    -o $path/bams/$muttype/$TUM_ID.neoantigens.chr$chr.$pos.$muttype.bam --samoutputformat BAM;
    count=$(samtools view $path/bams/$muttype/$TUM_ID.neoantigens.chr$chr.$pos.$muttype.bam | wc -l | awk '{print $1}');
    echo $chr $pos $ref $alt $gene chr$chr.$pos.$ref.$alt.$gene $count >> $path/counts/$TUM_ID.neoantigens.vcf.transcript.$muttype.gene.counts.tsv; # writing counts to file
done < $path/variants/$TUM_ID.neoantigens.$muttype.tsv
touch $path/logs/$TUM_ID.$muttype.extract_reads.success;
#<(tail -n "+$linenumber" $path/${TUM_ID}.neoantigens.${CLASS}.filtered.full.vcf.tsv)
