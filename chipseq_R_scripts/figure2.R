### This script is used to generate figures illustrating results from peak annotation using HOMER

## Load libraries
library(ggplot2)

## Set color palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
cols <- gg_color_hue(8)
cols[7] <- "#CC0033" #change to red

## SHSY5Y
allfiles <- scan(file="SHSY5Y-annotate-regions/filenames.txt", what=character())
listValues <- list()

for (n in 1:length(allfiles)){
  # Load in file
  filename <- allfiles[n]
  file <- read.delim(paste("SHSY5Y-annotate-regions/", filename, sep=""), header=T, sep="\t")
  
  # Get annotations
  annot <- as.character(file$Annotation)
  # removeText <- function(x){gsub( " *\\(.*?\\) *", "", x)}
  removeText <- function(x){strsplit(x, split=" \\(")[[1]][1]}
  file$Annotation <- factor(sapply(annot, removeText, USE.NAMES=F))
  
  # Save to list
  values <- (summary(file$Annotation)/nrow(file)) * 100
  label <- strsplit(filename, split="-")[[1]][1:2]
  label <- paste(label, collapse="_")
  listValues[[n]] <- cbind(names(values), rep(label, length(values)), values)
}


listValues[[2]] <- rbind(c("3' UTR", "SHSIR_ins", "0"), listValues[[2]])
df <- data.frame(do.call(rbind, listValues))
df$values <- as.numeric(as.character(df$values))
names(df) <- c("annotation", "sample", "values")

ggplot(df, aes(x=annotation, y=values)) +
  theme_bw(base_size=10) +
  theme(panel.grid.major = element_line(size = .5, color = "grey"),
        axis.text.x = element_text(angle=45, hjust=1)) +
  geom_bar(aes(fill=annotation),  position="dodge", stat="identity") +
  ylab("% of total regions") + xlab("") +
  theme(axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) +
  facet_grid(. ~ sample, scales="free") + 
  scale_fill_manual(values = cols) + 
  ylim(c(0,80))


## HepG2

allfiles <- scan(file="HepG2-annotate-regions/filenames.txt", what=character())
listValues <- list()

for (n in 1:length(allfiles)){
  # Load in file
  filename <- allfiles[n]
  file <- read.delim(paste("HepG2-annotate-regions/", filename, sep=""), header=T, sep="\t")
  
  # Get annotations
  annot <- as.character(file$Annotation)
  # removeText <- function(x){gsub( " *\\(.*?\\) *", "", x)}
  removeText <- function(x){strsplit(x, split=" \\(")[[1]][1]}
  file$Annotation <- factor(sapply(annot, removeText, USE.NAMES=F))
  
  # Save to list
  values <- (summary(file$Annotation)/nrow(file)) * 100
  label <- strsplit(filename, split="-")[[1]][1:2]
  label <- paste(label, collapse="_")
  listValues[[n]] <- cbind(names(values), rep(label, length(values)), values)
}

df <- data.frame(do.call(rbind, listValues))
df$values <- as.numeric(as.character(df$values))
names(df) <- c("annotation", "sample", "values")

ggplot(df, aes(x=annotation, y=values)) +
  theme_bw(base_size=10) +
  theme(panel.grid.major = element_line(size = .5, color = "grey"),
        axis.text.x = element_text(angle=45, hjust=1)) +
  geom_bar(aes(fill=annotation),  position="dodge", stat="identity") +
  ylab("% of total regions") + xlab("") +
  theme(axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) +
  facet_grid(. ~ sample, scales="free") +
  scale_fill_manual(values = cols) +
  ylim(c(0,80))


