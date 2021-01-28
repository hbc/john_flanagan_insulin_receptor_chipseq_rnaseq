module load gcc/6.2.0  python/2.7.12
module load beta/1.0.7

cd /n/data1/cores/bcbio/PIs/john_flanagan/flanagan-chipseq2rna/manuscript_revisions/BETA

#BETA plus -p beta-HepG2IR-ins-idr  -e res_treatment_ins_vs_con_noNA.bsf  -k BSF -g hg19 --gs /n/groups/bcbio/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa -o 100kb_genes_0.001/ -n BETA_IR_siggenes -d 100000 --df 0.001 --info 1,2,3 --gname2 -c 1

#BETA plus -p beta-ENCFF632WJV_H3K4me3.bed  -e res_treatment_ins_vs_con_noNA.bsf  -k BSF -g hg19 --gs /n/groups/bcbio/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa -o H3K4me3_ENCFF632WJV_100kb_genes_0.001/ -n BETA_H3K4me3_siggenes -d 100000 --df 0.001 --info 1,2,3 --gname2 -c 1

# Adding a parameter to change the distance we are searching

# IR
BETA plus -p beta-HepG2IR-ins-idr.5c.bed  -e res_treatment_ins_vs_con_noNA.bsf  -k BSF -g hg19 --gs /n/groups/bcbio/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa -o IR_1kb_genes_0.001/ -n BETA_IR_siggenes -d 1000 --df 0.001 --info 1,2,3 --gname2 -c 1

# Pol2
BETA plus -p beta-HepG2Pol2-ins-idr.5c.bed  -e res_treatment_ins_vs_con_noNA.bsf  -k BSF -g hg19 --gs /n/groups/bcbio/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa -o Pol2_1kb_genes_0.001/ -n BETA_Pol2_siggenes -d 1000 --df 0.001 --info 1,2,3 --gname2 -c 1

# Melissa's H3K4me3 sample
BETA plus -p beta-ENCFF632WJV_H3K4me3.bed.5c.bed  -e res_treatment_ins_vs_con_noNA.bsf  -k BSF -g hg19 --gs /n/groups/bcbio/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa -o H3K4me3_ENCFF632WJV_1kb_genes_0.001/ -n BETA_H3K4me3_siggenes -d 1000 --df 0.001 --info 1,2,3 --gname2 -c 1

# HCFC1
BETA plus -p beta-HCFC1-encode-idr.5c.bed  -e res_treatment_ins_vs_con_noNA.bsf  -k BSF -g hg19 --gs /n/groups/bcbio/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa -o HCFC1_1kb_genes_0.001/ -n BETA_HCFC1_siggenes -d 1000 --df 0.001 --info 1,2,3 --gname2 -c 1
