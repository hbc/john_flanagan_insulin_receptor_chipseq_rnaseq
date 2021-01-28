#! /bin/bash

module load  gcc/6.2.0  python/2.7.12
module load deeptools/2.5.3

for file in /n/data1/cores/bcbio/flanagan-chipseq2rna/chipseq-analysis-HepG2/bowtie1/*IR*.bam
do
base=`basename $file _R1.fastq_bowtie_sorted.bam`
bamCoverage -b $file -p 6 -bs 20 -o /n/data1/cores/bcbio/flanagan-chipseq2rna/chipseq-analysis-HepG2/bigWig/${base}.bw
done


