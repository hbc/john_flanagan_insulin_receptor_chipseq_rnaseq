## DiffBind

library(DiffBind)

# Read in samplesheet
samples <- read.delim("DiffBind_samplesheet.txt", header=T, sep="\t")

# Create dbObj
dbObj <- dba(sampleSheet=samples)

# Create affinity binding matrix
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)

# use dbOj to plot PCA
dba.plotPCA(dbObj,  attributes=DBA_CONDITION, label=DBA_ID)

# Plot heatmap
plot(dbObj)

# Establishing a contrast
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION, minMembers = 2)

# Performing differential analysis
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
dba.show(dbObj, bContrasts=T)	


# Overlap between methods
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)

# MA plot
dba.plotMA(dbObj, method=DBA_DESEQ, bNormalized = TRUE, bUsePval = FALSE, th=0.05)
dba.plotMA(dbObj, method=DBA_DESEQ, bNormalized = TRUE, bUsePval = FALSE, th=0.05, bXY = T)

# Write to file
# DESeq2

comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ, contrast = 1, th=1)

out <- as.data.frame(comp1.deseq)
write.table(out, file="diffbind/HepG2_ins_vs_con_deseq.txt", sep="\t", quote=F, col.names = NA)

out <- as.data.frame(comp1.edgeR)
write.table(out, file="diffbind/HepG2_ins_vs_con_edger.txt", sep="\t", quote=F, col.names = NA)


