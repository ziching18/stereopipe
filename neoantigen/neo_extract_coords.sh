TUM_ID=$1;
path=$(echo /stereoseq/all_samples/neoantigens/$TUM_ID);

for file in $(ls $path/bams/*.bam); do
    touch $path/coords/$TUM_ID.neoantigens.filtered.coords.txt;
    samtools view $file | awk '{print $1,$(NF-1),$NF}' >> $path/coords/$TUM_ID.neoantigens.filtered.coords.txt;
done