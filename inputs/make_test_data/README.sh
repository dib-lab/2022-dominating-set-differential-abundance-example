cd /group/ctbrowngrp2/tereiter/github/2020-ibd/outputs/sgc_genome_queries

for sample in PSM7J199 PSM7J1BJ MSM9VZMM MSM9VZL5 p8775mo1 p8816mo1 p9220mo1 p9223mo1
do
zcat ${sample}_k31_r1_search_oh0/GCF_008121495.1*reads* ${sample}_k31_r1_search_oh0/GCF_002234575.2*reads* > ~/github/2022-dominating-set-differential-abundance-example/inputs/make_test_data/${sample}.fq
done

cd ~/github/2022-dominating-set-differential-abundance-example/inputs/make_test_data/ 
gzip *fq

conda activate /home/tereiter/github/2022-dominating-set-differential-abundance-example/.snakemake/conda/dbc0297d0123f0227f129a9e4555a9e7

# repair and split reads into R1 and R2 with bbmap
for infile in *fq.gz
do
bn=$(basename $infile .fq.gz)
repair.sh in1=${infile} out1=${bn}_R1.fq.gz out2=${bn}_R2.fq.gz outs=${bn}_singleton.fq.gz repair
done
