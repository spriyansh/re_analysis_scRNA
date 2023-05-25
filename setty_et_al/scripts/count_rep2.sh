cellranger count \
	--id=rep2_p2 \
	--fastqs=../input/ \
	--sample=P2 \
	--transcriptome=../output/cr_index_hs/ \
	--localmem=62 \
	--localcores=8

# Move
mv rep2_p2 ../output/
