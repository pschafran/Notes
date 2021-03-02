#Use DADA to identify Isoetes PacBio LFY amplicons
# https://benjjneb.github.io/dada2/tutorial.html
# https://benjjneb.github.io/LRASManuscript/LRASms_fecal.html
## PacBio DADA Tutorial
library(dada2);packageVersion("dada2")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")
library(reshape2); packageVersion("reshape2")
library(gridExtra); packageVersion("gridExtra")
library(phyloseq); packageVersion("phyloseq")
library(DECIPHER)

setwd("~/Desktop/DADA/")
### After reads split and properly oriented continue here

path1 <- "~/Desktop/DADA/seqs/reoriented_seqs"
path.out <- "~/Desktop/DADA/Figures/"
path.rds <- "~/Desktop/DADA/RDS/"
theme_set(theme_bw())

fns1 <- list.files(path1, pattern="fastq", full.names=TRUE)
lens.fn <- lapply(fns1, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)

# Filter reads by length, quality and expected errors
filts1 <- file.path("filtered", basename(fns1))
track1 <- filterAndTrim(fns1, filts1, minQ=3, minLen=1000, maxLen=1200, trimLeft = 16, trimRight = 16, maxN=0, rm.phix=FALSE, maxEE=5)
track1

filtered.path <- "~/Desktop/DADA/filtered"
flts1 <- fns1 <- list.files(filtered.path, pattern="fastq", full.names=TRUE)
lens.fn <- lapply(flts1, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)

drp1 <- derepFastq(flts1, verbose=TRUE)
err1 <- learnErrors(drp1, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)
saveRDS(err1, file.path(path.rds, "PacBio_err1.rds"))
plotErrors(err1)
dd2 <- dada(drp1, err=err1, BAND_SIZE=32, multithread=TRUE, pool="pseudo")
#saveRDS(dd2, file.path(path.rds, "PacBio_dd2.rds"))
cbind(filtered=track1[,2], denoised=sapply(dd2, function(x) sum(x$denoised)))
st2 <- makeSequenceTable(dd2); dim(st2)
bim2 <- isBimeraDenovo(st2, minFoldParentOverAbundance=3.5, multithread=TRUE)
table(bim2)
sum(st2[,bim2])/sum(st2)
View(dd2)
for (sample in names(dd2)) {
  seqList <- getSequences(dd2[[sample]], collapse = TRUE)
  i <- 1
  while (i <= length(seqList)){
    seqName <- as.character(paste(">", sample, "_ASV", i, sep = ""))
    write(seqName, file="JSS1_DADA_ASVs.fasta", append = TRUE)
    write(seqList[[i]], file="JSS1_DADA_ASVs.fasta", append = TRUE)
    i <- i+1
      }
}


for (sample in names(dd2)) {
  seqList <- getSequences(dd2[[sample]], collapse = TRUE)
  write(sample, file="JSS1_DADA_ASVs.tsv", append = TRUE)
  write(length(seqList), file="JSS1_DADA_ASVs.tsv", append = TRUE)
}

