cellranger count \
	--id=rep1_p1 \
	--fastqs=../input/ \
	--sample=P1 \
	--transcriptome=../output/cr_index_hs/ \
	--localmem=62 \
	--localcores=8

# Move
mv rep1_p1 ../output/
