# Run All samples
cellranger count \
        --id=P3Comb \
        --fastqs=../input/ \
        --sample=P3Comb1,P3Comb2,P3Comb3,P3Comb4,P3Comb5,P3Comb6,P3Comb7,P3Comb8 \
        --transcriptome=../input2/cr_index_hs/ \
        --localmem=48 \
        --localcores=6

# Move
mv P3Comb ../output/