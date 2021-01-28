## DiffBind
library(DiffBind)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(tidyverse)


# Read in samplesheet
samples <- read.delim("DiffBind_samplesheet.txt", header=T, sep="\t")

# Create dbObj
dbObj <- dba(sampleSheet=samples)

# Create affinity binding matrix
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)

# Create several count matrices
dbObj_counts <- dba.count(dbObj, peaks=NULL,score=DBA_SCORE_READS_MINUS)
dbObj_full <- dba.count(dbObj, peaks=NULL,score=DBA_SCORE_TMM_MINUS_FULL)
dbObj_eff <- dba.count(dbObj, peaks=NULL,score=DBA_SCORE_TMM_MINUS_EFFECTIVE)
dbObj_chip_TMM <- dba.count(dbObj, peaks=NULL,score=DBA_SCORE_TMM_READS_FULL)

# Output a dataframe
counts <- dba.peakset(dbObj_counts, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
tmm_full <- dba.peakset(dbObj_full, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
tmm_eff <- dba.peakset(dbObj_eff, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
tmm_chip <- dba.peakset(dbObj_chip_TMM, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)

# Create a GenomicRanges object
gr_full <- GRanges(seqnames = tmm_full[,1], strand = rep("*", nrow(tmm_full)),
              ranges = IRanges(start = tmm_full[,2], end = tmm_full[,3]), 
                               mcols=tmm_full[,4:7])
gr_eff <- GRanges(seqnames = tmm_eff[,1], strand = rep("*", nrow(tmm_eff)),
                   ranges = IRanges(start = tmm_eff[,2], end = tmm_eff[,3]), 
                   mcols=tmm_eff[,4:7])

gr_chip <- GRanges(seqnames = tmm_chip[,1], strand = rep("*", nrow(tmm_chip)),
                  ranges = IRanges(start = tmm_chip[,2], end = tmm_chip[,3]), 
                  mcols=tmm_chip[,4:7])
# Set database
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Annotate those peaks
tmm_full_annot <- annotatePeak(gr_full, tssRegion = c(-2000, 2000), TxDb = txdb)
tmm_eff_annot <- annotatePeak(gr_eff, tssRegion = c(-2000, 2000), TxDb = txdb)
tmm_chip_annot <- annotatePeak(gr_chip, tssRegion = c(-2000, 2000), TxDb = txdb)

# Need to go from Entrez to Ensembl/Gene symbols
# Get annotations
require(AnnotationHub)
ah <- AnnotationHub()
orgdb <- query(ah, c("Homo sapiens", "EnsDb"))[["AH60977"]]

# See what columns we have and get what we need
keytypes(orgdb)
egid <- as.character(keys(orgdb, "ENTREZID"))
ensAnnot <- select(orgdb, egid, c("SYMBOL", "GENEID"), "ENTREZID")
ensAnnot$ENTREZID <- as.character(ensAnnot$ENTREZID)
              

# Add columns
tmm_full_annot <- tmm_full_annot@anno %>% data.frame()
tmm_eff_annot <- tmm_eff_annot@anno %>% data.frame()
tmm_chip_annot <- tmm_chip_annot@anno %>% data.frame()

norm_counts_full <- tmm_full_annot %>% 
  data.frame() %>% 
  left_join(ensAnnot, by=c("geneId"="ENTREZID")) %>%
  select(seqnames, start, end, mcols.expt1_con, mcols.expt2_con, 
         mcols.expt1_ins, mcols.expt2_ins, SYMBOL, GENEID)

norm_counts_eff <- tmm_eff_annot %>% 
  data.frame() %>% 
  left_join(ensAnnot, by=c("geneId"="ENTREZID")) %>%
  select(seqnames, start, end, mcols.expt1_con, mcols.expt2_con, 
         mcols.expt1_ins, mcols.expt2_ins, SYMBOL, GENEID)

norm_counts_chip <- tmm_chip_annot %>% 
  data.frame() %>% 
  left_join(ensAnnot, by=c("geneId"="ENTREZID")) %>%
  dplyr::select(seqnames, start, end, mcols.expt1_con, mcols.expt2_con, 
         mcols.expt1_ins, mcols.expt2_ins, SYMBOL, GENEID)

# Subset these matrices to the up and down regulated genes

upgenes <- read.csv("~/Dropbox (HBC)/HBC Team Folder (1)/Consults/john_flanagan/flanagan-rnaseq/results/differential_expression/res_treatment_ins_vs_con_deg_lfc_up.csv.gz")

downgenes <- read.csv("~/Dropbox (HBC)/HBC Team Folder (1)/Consults/john_flanagan/flanagan-rnaseq/results/differential_expression/res_treatment_ins_vs_con_deg_lfc_down.csv.gz")

length(which(norm_counts_chip$GENEID %in% upgenes$ensgene))
length(which(norm_counts_chip$GENEID %in% downgenes$ensgene))

# Write to file
out_full <- norm_counts_full[which(norm_counts_full$GENEID %in% upgenes$ensgene),]
write.csv(out_full, "files_for_boxplots/tmm_full_norm_upgenes.csv", row.names=FALSE, quote=F)

out_full_down <- norm_counts_full[which(norm_counts_full$GENEID %in% downgenes$ensgene),]
write.csv(out_full_down, "files_for_boxplots/tmm_full_norm_downgenes.csv", row.names=FALSE, quote=F)

out_eff_up <- norm_counts_eff[which(norm_counts_eff$GENEID %in% upgenes$ensgene),]
write.csv(out_eff_up, "files_for_boxplots/tmm_eff_norm_upgenes.csv", row.names=FALSE, quote=F)

out_eff_down <- norm_counts_eff[which(norm_counts_eff$GENEID %in% downgenes$ensgene),]
write.csv(out_eff_down, "files_for_boxplots/tmm_eff_norm_downgenes.csv", row.names=FALSE, quote=F)

out_chip_up <- norm_counts_chip[which(norm_counts_chip$GENEID %in% upgenes$ensgene),]
write.csv(out_chip_up, "files_for_boxplots/tmm_chip_norm_upgenes.csv", row.names=FALSE, quote=F)

out_chip_down <- norm_counts_chip[which(norm_counts_chip$GENEID %in% downgenes$ensgene),]
write.csv(out_chip_down, "files_for_boxplots/tmm_chip_norm_downgenes.csv", row.names=FALSE, quote=F)
