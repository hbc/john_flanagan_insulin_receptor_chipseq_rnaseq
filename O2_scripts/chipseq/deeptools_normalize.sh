#! /bin/bash

module load  gcc/6.2.0  python/2.7.12
module load deeptools/2.5.3

cd /n/data1/cores/bcbio/PIs/john_flanagan/flanagan-chipseq2rna/manuscript_revisions/deepTools 


# Pol2 control

bamCompare -b1 ../../chipseq-analysis-HepG2/bowtie1/20150624-7-HepG2-con-2-Pol2-MH2026-1_S7_R1.fastq_bowtie_sorted.bam \
-b2 ../../chipseq-analysis-HepG2/bowtie1/20150624-11-HepG2-con-2-input-MH2026-1_S11_R1.fastq_bowtie_sorted.bam \
-o normalized_to_input_bw/HepG2-con-2-Pol2.bw \
--binSize 50 \
--normalizeTo1x 130000000 \
--smoothLength 60 \
--extendReads 150 \
--centerReads \
-p 6 2> logs/HepG2-con-2-Pol2_bamCompare.log


# Pol2 ins

bamCompare -b1 ../../chipseq-analysis-HepG2/bowtie1/20150624-8-HepG2-ins-2-Pol2-MH2026-1_S8_R1.fastq_bowtie_sorted.bam \
-b2 ../../chipseq-analysis-HepG2/bowtie1/20150624-12-HepG2-ins-2-input-MH2026-1_S12_R1.fastq_bowtie_sorted.bam \
-o normalized_to_input_bw/HepG2-ins-2-Pol2.bw \
--binSize 50 \
--normalizeTo1x 130000000 \
--smoothLength 60 \
--extendReads 150 \
--centerReads \
-p 6 2> logs/HepG2-ins-2-Pol2_bamCompare.log


# IR control

bamCompare -b1 ../../chipseq-analysis-HepG2/bowtie1/20150624-3-HepG2-con-1-IR-MH2026-1_S3_R1.fastq_bowtie_sorted.bam \
-b2 ../../chipseq-analysis-HepG2/bowtie1/20150624-5-HepG2-con-1-input-MH2026-1_S5_R1.fastq_bowtie_sorted.bam \
-o normalized_to_input_bw/HepG2-con-1-IR.bw \
--binSize 50 \
--normalizeTo1x 130000000 \
--smoothLength 60 \
--extendReads 150 \
--centerReads \
-p 6 2> logs/HepG2-con-1-IR_bamCompare.log


# IR ins

bamCompare -b1 ../../chipseq-analysis-HepG2/bowtie1/20150624-4-HepG2-ins-1-IR-MH2026-1_S4_R1.fastq_bowtie_sorted.bam \
-b2 ../../chipseq-analysis-HepG2/bowtie1/20150624-6-HepG2-ins-1-input-MH2026-1_S6_R1.fastq_bowtie_sorted.bam \
-o normalized_to_input_bw/HepG2-ins-1-IR.bw \
--binSize 50 \
--normalizeTo1x 130000000 \
--smoothLength 60 \
--extendReads 150 \
--centerReads \
-p 6 2> logs/HepG2-ins-1-IR_bamCompare.log

# HCFC1

bamCompare -b1 ../../encode/HCFC1/bams/ENCFF002EDN.sorted.bam \
-b2 ../../encode/HCFC1/bams/ENCFF002ECU.sorted.bam \
-o normalized_to_input_bw/HCFC1_ENCFF002EDN.bw \
--binSize 50 \
--normalizeTo1x 130000000 \
--smoothLength 60 \
--extendReads 150 \
--centerReads \
-p 6 2> logs/HCFC1_ENCFF002EDN_bamCompare.log



bamCompare -b1 ../../encode/HCFC1/bams/ENCFF002EDS.sorted.bam \
-b2 ../../encode/HCFC1/bams/ENCFF002ECQ.sorted.bam \
-o normalized_to_input_bw/HCFC1_ENCFF002EDS.bw \
--binSize 50 \
--normalizeTo1x 130000000 \
--smoothLength 60 \
--extendReads 150 \
--centerReads \
-p 6 2> logs/HCFC1_ENCFF002EDS_bamCompare.log


