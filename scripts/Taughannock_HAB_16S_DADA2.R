setwd("~/Desktop/HAB/16S_barcoding/")
library(dada2)
library(ggplot2)
library(viridis)
library(RColorBrewer)
path <- "~/Desktop/HAB/16S_barcoding/"
list.files(path)

# Get sample files
fnFs <- sort(list.files(path, pattern="_R1_fastp.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_fastp.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# View sequence qualities
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), trimLeft = c(20, 20),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# Learn error profile of the data
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)


# Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# ASV inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merge sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Identify chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Output final ASV seqs
uniquesToFasta(seqtab.nochim, "HAB_16S_ASVs.fasta")

# Track Reads QC
getN <- function(x) sum(getUniques(x))
track <- cbind(out, getN(dadaFs), getN(dadaRs), getN(mergers), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Assign taxonomy
silva.taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE, minBoot = 50, tryRC = TRUE)
silva.taxa <- addSpecies(silva.taxa, "silva_species_assignment_v138.1.fa.gz", verbose = TRUE)
silva.taxa.print <- silva.taxa # Removing sequence rownames for display only
rownames(silva.taxa.print) <- NULL
View(silva.taxa.print)
View(silva.taxa)

rdp.taxa <- assignTaxonomy(seqtab.nochim, "rdp_train_set_18.fa.gz", multithread=TRUE, minBoot = 50, tryRC = TRUE)
rdp.taxa <- addSpecies(rdp.taxa, "rdp_species_assignment_18.fa.gz", verbose = TRUE)
rdp.taxa.print <-rdp.taxa
rownames(rdp.taxa.print) <- NULL
View(rdp.taxa.print) 


phylum_asvs <- data.frame(Dummy = 1, Count = c(478,18,7022,7,10,16,9507))
rownames(phylum_asvs) <- c("Bacteriodota","Bdellovibrionota","Cyanobacteria","Deinococcota","Firmicutes","Myxococcota","Proteobacteria")

phylum_asvs_plot <- ggplot(phylum_asvs) +
  geom_bar(aes_(y=phylum_asvs$Count, x= phylum_asvs$Dummy, fill=rownames(phylum_asvs)), position="fill", stat="identity", width = 0.5) +
  scale_fill_manual(values=rep(brewer.pal(12,"Paired"),times=40)) +
  labs(x="", y="Perc. annotated ASVs", fill = "")
phylum_asvs_plot + theme_bw()




class_asvs <- data.frame(Dummy = 1, Count = c(446,32,18,7022,7,10,16,5558,3949))
rownames(class_asvs) <- c("Bacteroidia","Kapabacteria","Oligoflexia","Cyanobacteria","Deinococci","Bacilli","Polyangia","Alphaproteobacteria","Gammaproteobacteria")

class_asvs_plot <- ggplot(class_asvs) +
  geom_bar(aes_(y=class_asvs$Count, x= class_asvs$Dummy, fill=rownames(class_asvs)), position="fill", stat="identity", width = 0.5) +
  scale_fill_manual(values=rep(brewer.pal(12,"Paired"),times=40)) +
  labs(x="", y="Perc. annotated ASVs", fill = "")
class_asvs_plot + theme_bw()




order_asvs <- data.frame(Dummy = 1, Count = c(15,155,271,2,32,18,7007,11,6,7,10,6,10,590,962,223,42,402,21,10,2933,3352,569,28))
rownames(order_asvs) <- c("Chitinophagales","Cytophagales","Flavobacteriales","Spingobacteriales","Kapabacteriales","Silvanigrellales","Cyanobacteriales","Pseudanabaenales","Synechococcales","Deinococcales","Exiguobacterales","Blfdi19","Polyangiales","Caulobacterales","Elsterales","Rhizobiales","Rhodobacterales","Rhodospirillales","SAR11 clade","Rickettsiales","Sphingomonadales","Burkholderiales","Enterobacteriales","Pseudomonales")

order_asvs_plot <- ggplot(order_asvs) +
  geom_bar(aes_(y=order_asvs$Count, x= order_asvs$Dummy, fill=rownames(order_asvs)), position="fill", stat="identity", width = 0.5) +
  scale_fill_manual(values=rep(brewer.pal(12,"Paired"),times=40)) +
  labs(x="", y="Perc. annotated ASVs", fill = "")
order_asvs_plot + theme_bw()





family_asvs <- data.frame(Dummy = 1, Count = c(15,25, 128, 271, 2, 18, 2, 6989, 11, 6, 7, 10, 10, 87, 503, 962, 67, 47, 105, 42, 402, 5, 2933, 9, 19, 20, 3293, 11, 373, 196, 8, 20))
rownames(family_asvs) <- c("Chitinophagaceae", "Bernardetiaceae","Spirosomaceae","Flabocteriaceae","LiUU-11-161","Silvanigrellaceae","Microcystaceae","Nostocaceae","Pseudanabaenaceae","Cyanobiaceae","Deinococcaceae","Exiguobacteraceae","Polyangiaceae","Caulobacteraceae","Hyphomonadaceae","Elsteraceae","Beijerinckiaceae","Rhizobiaceae","Rhizobiales Incertae Sedis","Rhodobacteraceae","Rhodospirillaceae","Rickettsiaceae","Sphingomonadaceae","Chitinibacteriaceae","Chitinimonadaceae","Chromobacteriaceae","Comamonadaceae","Rhodocyclaceae","Aeromonadaceae","Altermonadaceae","Cellvibrionaceae","Pseudomonadaceae")

family_asvs_plot <- ggplot(family_asvs) +
  geom_bar(aes_(y=family_asvs$Count, x= family_asvs$Dummy, fill=rownames(family_asvs)), position="fill", stat="identity", width = 0.5) +
  scale_fill_manual(values=rep(brewer.pal(12,"Paired"),times=40)) +
  labs(x="", y="Perc. annotated ASVs", fill = "")
family_asvs_plot + theme_bw()






genus_asvs <- data.frame(Dummy = 1, Count= c(6,9,25,122,6,271,2,6989,11,6,7,10,10,27,3,57,503,382,4,67,4,47,105,402,5,2933,9,19,20,269,261,79,186,1952,373,196,8,20))
rownames(genus_asvs) <- c("Ferruginibacter","Sediminibacterium","Bernardetia","Flectobacillus","Lacihabitans","Flavobacterium","Microcystis","Aphanizomenon","Pseudanabaena","Cyanobium","Deinococcus","Exiguobacterium","Pajaroellobacter","Brevundimonas","Caulobacter","Phenylobacterium","UKL13-1","Elstera","Lacibacterium","Bosea","FukuN57","Rhizobium s.l.", "Phreatobacter","Roseospirillum","Candidatus Megaira","Sphingorhabus","Deefgea","Chitinimonas","Vogesella","Hydrogenophaga","Inhella","Leptothrix","Pelomonas","Rhodoferax","Aeromonas","Rheinheimera","Cellvibrio","Pseudomonas")

genus_asvs_plot <- ggplot(genus_asvs) +
  geom_bar(aes_(y=genus_asvs$Count, x= genus_asvs$Dummy, fill=rownames(genus_asvs)), position="fill", stat="identity", width = 0.5) +
  scale_fill_manual(values=rep(brewer.pal(12,"Paired"),times=40)) +
  labs(x="", y="Perc. annotated ASVs", fill = "")
genus_asvs_plot + theme_bw()



