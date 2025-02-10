TUM_ID=$1
muttype=$2;
path=$(echo /stereoseq/all_samples/mutations/$TUM_ID);

for file in $(ls $path/bams/$muttype/*.$muttype.bam); do
    touch $path/coords/$TUM_ID.somatic.$muttype.coords.txt;
    samtools view $file | awk '{print $1,$(NF-1),$NF}' >> $path/coords/$TUM_ID.somatic.$muttype.coords.txt;
done
