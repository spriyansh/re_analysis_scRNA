# Run
cellranger mkref \
	--genome=cr_index_hs \
	--fasta=../input/GRCh38.primary_assembly.genome.fa \
	--genes=../output/gencode.v43.annotatio.filtered.gtf

# Move
mv cr_index_hs ../output/
