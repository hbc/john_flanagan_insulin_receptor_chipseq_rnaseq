
#!/bin/sh

date

module load seq/macs/2.1.0

#File info
#fastqFile=20150319-1-HepG2input-MH1746_S1_R1.fastq.gz
fastqFile=$1

#Directories
fastqDir=/n/data1/cores/bcbio/flanagan-chipseq2rna/HepG2_chipseq
baseDir=/n/data1/cores/bcbio/flanagan-chipseq2rna/chipseq-analysis/bowtie1


echo $fastqFile

#Trimming
/groups/bcbio/bcbio/anaconda/bin/cutadapt --times=2 --quality-base=33 --quality-cutoff=5 --format=fastq --adapter=AGATCGGAAGAG --adapter=CTCTTCCGATCT  --minimum-length=25  -o $baseDir/${fastqFile}.trimmed.gz $fastqDir/$fastqFile

gzip -d $baseDir/${fastqFile}.trimmed.gz

#Mapping trimmed reads
bowtie -q -S  -m 1 -p 8 /groups/shared_databases/bowtie_indexes/hg19 $baseDir/${fastqFile}.trimmed $baseDir/${fastqFile}_bowtie.sam

# Converting to bam
/opt/bcbio/local/bin/samtools view -@ 8 -h -S -b $baseDir/${fastqFile}_bowtie.sam -o $baseDir/${fastqFile}_bowtie.bam

# Sorting
 /opt/bcbio/local/bin/sambamba sort -t 8  -o $baseDir/${fastqFile}_bowtie_sorted.bam $baseDir/${fastqFile}_bowtie.bam


#Indexing
samtools index $baseDir/${fastqFile}_bowtie_sorted.bam

#Clean up
#rm $baseDir/${fastqFile}.trimmed
rm $baseDir/${fastqFile}_bowtie.sam
rm $baseDir/${fastqFile}_bowtie.bam

#Peak calling
#macs2 callpeak -t $baseDir/${fastqFile}_bowtie_sorted.bam  -f BAM -g hs -n ${fastqFile}_bowtie 2> $baseDir/${fastqFile}_bowtie_macs2.log

date

