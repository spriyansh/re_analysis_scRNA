cellranger count \
	--id=rep3_p3 \
	--fastqs=../input/ \
	--sample=P3 \
	--transcriptome=../output/cr_index_hs/ \
	--localmem=62 \
	--localcores=8

# Move
mv rep3_p3 ../output/
