## Generating heatmaps around the TSS


# Loading packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)

# Load files
HepG2 <- list("HOMER//heatmaps_manuscript//HepG2IR-IP-ins_peaks.bed", "HOMER/heatmaps_manuscript/HepG2Pol2-IP-ins_peaks.bed")
names(HepG2) <- c("HepG2IR", "HepG2Pol2")

peakIR <- readPeakFile(HepG2[[1]])
peakPol2 <- readPeakFile(HepG2[[2]])

# Heatmap IR
promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)
tagMatrix <- getTagMatrix(peakIR, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-2000, 2000), color="red")

# Heatmap Pol2
#promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)
tagMatrix <- getTagMatrix(peakPol2, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-2000, 2000), color="red")

# Load files
SHS<- list("HOMER//heatmaps_manuscript//SHSY5YIR-ins1_peaks.bed", "HOMER/heatmaps_manuscript/SHSY5YPol2-ins1_peaks.bed")
names(SHS) <- c("SHSIR", "SHSPol2")

peakIR <- readPeakFile(SHS[[1]])
peakPol2 <- readPeakFile(SHS[[2]])

# Heatmap IR
promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)
tagMatrix <- getTagMatrix(peakIR, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-2000, 2000), color="red")

# Heatmap Pol2
#promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)
tagMatrix <- getTagMatrix(peakPol2, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-2000, 2000), color="red")
