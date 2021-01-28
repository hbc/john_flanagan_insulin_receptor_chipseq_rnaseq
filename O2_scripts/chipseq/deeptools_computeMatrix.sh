#!pG2 /bin/bash

module load  gcc/6.2.0  python/2.7.12
module load deeptools/2.5.3

cd /n/data1/cores/bcbio/PIs/john_flanagan/flanagan-chipseq2rna/manuscript_revisions/deepTools 

#computeMatrix reference-point -S ../../chipseq-analysis-HepG2/bigWig/20150624-1-HepG2-con-1-Pol2-MH2026-1_S1.bw ../../chipseq-analysis-HepG2/bigWig/20150624-2-HepG2-ins-1-Pol2-MH2026-1_S2.bw  --referencePoint TSS -R upgene.bed  -b 2000 -a 4000 --skipZeros -p 6 -o Pol2_ins_vs_con_replicate1_upgenes.mat.gz

#computeMatrix reference-point -S ../../chipseq-analysis-HepG2/bigWig/20150624-1-HepG2-con-1-Pol2-MH2026-1_S1.bw ../../chipseq-analysis-HepG2/bigWig/20150624-2-HepG2-ins-1-Pol2-MH2026-1_S2.bw  --referencePoint TSS -R downgene.bed -b 2000 -a 4000 --skipZeros -p 6 -o Pol2_ins_vs_con_replicate1_downgenes.mat.gz

#computeMatrix reference-point -S ../../chipseq-analysis-HepG2/bigWig/20150624-7-HepG2-con-2-Pol2-MH2026-1_S7.bw ../../chipseq-analysis-HepG2/bigWig/20150624-8-HepG2-ins-2-Pol2-MH2026-1_S8.bw  --referencePoint TSS -R upgene.bed  -b 2000 -a 4000 --skipZeros -p 6 -o Pol2_ins_vs_con_replicate2_upgenes.mat.gz

#computeMatrix reference-point -S ../../chipseq-analysis-HepG2/bigWig/20150624-7-HepG2-con-2-Pol2-MH2026-1_S7.bw ../../chipseq-analysis-HepG2/bigWig/20150624-8-HepG2-ins-2-Pol2-MH2026-1_S8.bw  --referencePoint TSS -R downgene.bed -b 2000 -a 4000 --skipZeros -p 6 -o Pol2_ins_vs_con_replicate2_downgenes.mat.gz

#plotProfile -m Pol2_ins_vs_con_replicate1.mat.gz -out Pol2_ins_vs_con_rep1_profile.png --perGroup --plotTitle "" --samplesLabel "Control" "Insulin"

#plotProfile -m Pol2_ins_vs_con_replicate1_upgenes.mat.gz -out Pol2_ins_vs_con_rep1_upgenes_profile.png --perGroup --plotTitle "" --samplesLabel "Control" "Insulin" --regionsLabel "Up-regulated genes"

#plotProfile -m Pol2_ins_vs_con_replicate1_downgenes.mat.gz -out Pol2_ins_vs_con_rep1_downgenes_profile.png --perGroup --plotTitle "" --samplesLabel "Control" "Insulin" --regionsLabel "Down-regulated genes"

#computeMatrix reference-point -S normalized_to_input_bw/HepG2-con-2-Pol2.bw normalized_to_input_bw/HepG2-ins-2-Pol2.bw  --referencePoint TSS -R upgene.bed  -b 2000 -a 2000 --skipZeros -p 6 -o Pol2_ins_vs_con_replicate2_upgenes_norm.mat.gz

computeMatrix reference-point -S normalized_to_input_bw/HepG2-con-1-IR.bw normalized_to_input_bw/HepG2-ins-1-IR.bw  --referencePoint TSS -R upgene.bed  -b 2000 -a 2000 --skipZeros -p 6 -o IR_ins_vs_con_replicate1_upgenes_norm.mat.gz
