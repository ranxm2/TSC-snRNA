cellranger="/mnt/newhome/liangjie/software/cellranger-7.2.0/bin/cellranger"
ref="/mnt/newhome/liangjie/refFile/cellRanger/refdata-gex-GRCh38-2020-A"
fastqdir="/mnt/newhome/liangjie/WeiboNiu_Data/30-921109009"

for sample in CTRL PSZ-6 TSC-edge TSC-tuber
do
	echo $sample;
	out=$sample
	$cellranger count \
		--sample $sample \
		--id $out \
		--transcriptome $ref \
		--fastqs $fastqdir \
		--localcores 48 \
		--localmem 256 \
		--nosecondary \
		--disable-ui \
		--chemistry 'ARC-v1'
done

echo Job Done! 
