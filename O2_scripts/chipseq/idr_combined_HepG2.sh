#!/bin/sh

date

inputFile1=$1
treatFile1=$2
inputFile2=$3
treatFile2=$4
EXPT=$5

NAME1=${treatFile1##*/}
NAME1=${NAME1%-MH2026*}

NAME2=${treatFile2%-MH1746*}

#Directories
expt1Dir=/n/data1/cores/bcbio/flanagan-chipseq2rna/expt1-HepG2
macsDir=/n/data1/cores/bcbio/flanagan-chipseq2rna/idr_pipeline/idr_HepG2/macs
outputDir=/n/data1/cores/bcbio/flanagan-chipseq2rna/idr_pipeline/idr_HepG2/pooledPseudoReps
baseDir=/n/data1/cores/bcbio/flanagan-chipseq2rna/idr_pipeline/

#Peak calling

echo "Calling peaks for "$treatFile2
macs2 callpeak -t $expt1Dir/${treatFile2} -c $expt1Dir/${inputFile2} -f BAM -g hs -n $macsDir/${NAME2} -B -p 1e-3  2> $macsDir/${NAME2}_macs2.log

#Sort peak by -log10(p-value)
echo "Sorting peaks..."
sort -k8,8nr $macsDir/${NAME2}_peaks.narrowPeak | head -n 100000 > $macsDir/${NAME2}_sorted.narrowPeak

# Need Python3 local install
module unload dev/python/2.7.6

#Independent replicate IDR
echo "Running IDR on replicates..."
idr --samples ${baseDir}/macs/${NAME1}_sorted.narrowPeak ${macsDir}/${NAME2}_sorted.narrowPeak --input-file-type narrowPeak --output-file idr_HepG2/${EXPT}-idr --rank p.value --plot

#Merge treatment BAMS
echo "Merging BAM files for pseudoreplicates..."
samtools merge -u ${outputDir}/merged_bams/${NAME1}_${NAME2}_merged.bam /n/data1/cores/bcbio/flanagan-chipseq2rna/chipseq-analysis-HepG2/bowtie1/${treatFile1} $expt1Dir/${treatFile2}
samtools view -H ${outputDir}/merged_bams/${NAME1}_${NAME2}_merged.bam > ${EXPT}_header.sam

#Split merged treatments
nlines=$(samtools view ${outputDir}/merged_bams/${NAME1}_${NAME2}_merged.bam | wc -l ) # Number of reads in the BAM file
nlines=$(( (nlines + 1) / 2 )) # half that number
samtools view ${outputDir}/merged_bams/${NAME1}_${NAME2}_merged.bam | shuf - | split -d -l ${nlines} - "${outputDir}/${EXPT}" # This will shuffle the lines in the file and split it into two SAM files
cat ${EXPT}_header.sam ${outputDir}/${EXPT}00 | samtools view -bS - > ${outputDir}/${EXPT}00.bam
cat ${EXPT}_header.sam ${outputDir}/${EXPT}01 | samtools view -bS - > ${outputDir}/${EXPT}01.bam

#Merge input BAMS
echo "Merging input BAM files for pseudoreplicates..."
samtools merge -u ${outputDir}/merged_bams/${NAME1}input_${NAME2}input_merged.bam /n/data1/cores/bcbio/flanagan-chipseq2rna/chipseq-analysis-HepG2/bowtie1/${inputFile1} $expt1Dir/${inputFile2}

#Split merged treatment BAM
nlines=$(samtools view ${outputDir}/merged_bams/${NAME1}input_${NAME2}input_merged.bam | wc -l ) # Number of reads in the BAM file
nlines=$(( (nlines + 1) / 2 )) # half that number
samtools view ${outputDir}/merged_bams/${NAME1}input_${NAME2}input_merged.bam | shuf - | split -d -l ${nlines} - "${outputDir}/${EXPT}_input" # This will shuffle the lines in the file and split in two 
cat ${EXPT}_header.sam ${outputDir}/${EXPT}_input00 | samtools view -bS - > ${outputDir}/${EXPT}_input00.bam
cat ${EXPT}_header.sam ${outputDir}/${EXPT}_input01 | samtools view -bS - > ${outputDir}/${EXPT}_input01.bam

# Need Python2.7 for MACS2
module load dev/python/2.7.6

#Peak calling on pseudoreplicates
echo "Calling peaks for pseudoreplicate1 "
macs2 callpeak -t ${outputDir}/${EXPT}00.bam -c ${outputDir}/${EXPT}_input00.bam -f BAM -g hs -n $macsDir/${NAME1}_pr -B -p 1e-3  2> $macsDir/${NAME1}_pr_macs2.log

#echo "Calling peaks for pseudoreplicate2"
macs2 callpeak -t ${outputDir}/${EXPT}01.bam -c ${outputDir}/${EXPT}_input01.bam -f BAM -g hs -n $macsDir/${NAME2}_pr -B -p 1e-3  2> $macsDir/${NAME2}_pr_macs2.log

#Sort peak by -log10(p-value)
echo "Sorting peaks..."
sort -k8,8nr $macsDir/${NAME1}_pr_peaks.narrowPeak | head -n 100000 > $macsDir/${NAME1}_pr_sorted.narrowPeak
sort -k8,8nr $macsDir/${NAME2}_pr_peaks.narrowPeak | head -n 100000 > $macsDir/${NAME2}_pr_sorted.narrowPeak

# Need Python3 local install
module unload dev/python/2.7.6

#Independent replicate IDR
echo "Running IDR on pseudoreplicates..."
idr --samples ${macsDir}/${NAME1}_pr_sorted.narrowPeak ${macsDir}/${NAME2}_pr_sorted.narrowPeak --input-file-type narrowPeak --output-file idr_HepG2/${EXPT}_pseudorep-idr --rank p.value --plot


date

