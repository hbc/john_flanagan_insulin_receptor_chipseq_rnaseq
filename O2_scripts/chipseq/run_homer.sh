## First make tag directories using only the BED files
makeTagDirectory HCFC1_in_HepG2IR/ test_HCFC1_IR.bed -format bed

## Then use those tag directories to generate data for histograms (low mag)
 annotatePeaks.pl tss hg19 -size 4000 -hist 2 -d tagfiles_HCFC1_in_HepG2IR/ tagfiles_HepG2IR_in_HCFC1/  > HCFC1_IR.txt
