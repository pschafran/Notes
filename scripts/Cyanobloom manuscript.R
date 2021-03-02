library(ggplot2)
library(gridExtra)
library(viridis)
library(RColorBrewer)
setwd("~/Box/cyanobloom")


############# Taxonomic Comparisons b/w tools
{
genusDf <- read.table("~/Box/cyanobloom/classifier_superkingdom.table", header = TRUE, row.names = 1 , sep = "\t")
colnames(genusDf) <- c("Blast","Centrifuge", "Kaiju","Kraken","Agreement")
View(genusDf)
unique(genusDf$Blast_Genus)
genus_blast <- as.data.frame(table(genusDf$Blast, useNA = "always"), stringsAsFactors = F)
genus_centrifuge <- as.data.frame(table(genusDf$Centrifuge, useNA = "always"), stringsAsFactors = F)
genus_kaiju <- as.data.frame(table(genusDf$Kaiju, useNA = "always"), stringsAsFactors = F)
genus_kraken <- as.data.frame(table(genusDf$Kraken, useNA = "always"), stringsAsFactors = F)
genus_blast
genus_centrifuge
genus_kaiju
genus_kraken
mergedTableBC <- merge(genus_blast, genus_centrifuge, by = 1, all = TRUE)
mergedTableKK <- merge(genus_kaiju, genus_kraken, by = 1, all = TRUE)
mergedTable <- merge(mergedTableBC, mergedTableKK, by = 1, all = TRUE)
colnames(mergedTable) <- c("Superkingdom","Blast","Centrifuge", "Kaiju","Kraken")
View(mergedTable)

# start converting dataframe to format for ggplot
rowNames <- mergedTable[,1]
rowNamesList <- as.vector(rowNames)
rowNamesList[is.na(rowNamesList)] <- "Unclassified"
plotDf <- as.data.frame(mergedTable[,-1])
rownames(plotDf) <- rowNamesList
plotDf[is.na(plotDf)] <- 0
View(plotDf)

noUnclassifiedDf <- plotDf[!(row.names(plotDf)) %in% c("Unclassified", "root"), ]
blastGGdf <- as.data.frame(noUnclassifiedDf$Blast)
blastGGdf[,2] <- "Blast"
blastGGdf[,3] <- rownames(noUnclassifiedDf)
colnames(blastGGdf) <- c("count", "classifier", "genus")

centGGdf <- as.data.frame(noUnclassifiedDf$Centrifuge)
centGGdf[,2] <- "Centrifuge"
centGGdf[,3] <- rownames(noUnclassifiedDf)
colnames(centGGdf) <- c("count", "classifier", "genus")

kaijuGGdf <- as.data.frame(noUnclassifiedDf$Kaiju)
kaijuGGdf[,2] <- "Kaiju"
kaijuGGdf[,3] <- rownames(noUnclassifiedDf)
colnames(kaijuGGdf) <- c("count", "classifier", "genus")

krakenGGdf <- as.data.frame(noUnclassifiedDf$Kraken)
krakenGGdf[,2] <- "Kraken2"
krakenGGdf[,3] <- rownames(noUnclassifiedDf)
colnames(krakenGGdf) <- c("count", "classifier", "genus")

mergedGGBC <- rbind(blastGGdf, centGGdf)
mergedGGKK <- rbind(kaijuGGdf, krakenGGdf)
superkingdomMergedGGdf <- rbind(mergedGGBC, mergedGGKK)
#View(mergedGGdf)
}

superkingdomNum <- ggplot(data = superkingdomMergedGGdf) + 
  geom_bar(data = superkingdomMergedGGdf, aes(y = count, x = classifier, fill = genus), position="stack", stat="identity", width = 0.5) + 
  scale_fill_manual(values=rep(brewer.pal(12,"Paired"),times=5)) + 
  labs(x="", y="Num. of contigs", fill = "") +
  theme_classic(base_size = 25)
superkingdomNum

superkingdomProp <- ggplot(superkingdomMergedGGdf)+ 
  geom_bar(data = superkingdomMergedGGdf, aes(y = count, x = classifier, fill = genus), position="fill", stat="identity") + 
  scale_fill_manual(values=rep(brewer.pal(12,"Paired"),times=5)) + 
  labs(x="", y="Prop. annotated contigs", fill = "") +
  theme_classic(base_size = 25)
superkingdomProp

##### Phylum
{
  genusDf <- read.table("~/Box/cyanobloom/classifier_phylum.table", header = TRUE, row.names = 1 , sep = "\t")
  colnames(genusDf) <- c("Blast","Centrifuge", "Kaiju","Kraken","Agreement")
  genus_blast <- as.data.frame(table(genusDf$Blast, useNA = "always"), stringsAsFactors = F)
  genus_centrifuge <- as.data.frame(table(genusDf$Centrifuge, useNA = "always"), stringsAsFactors = F)
  genus_kaiju <- as.data.frame(table(genusDf$Kaiju, useNA = "always"), stringsAsFactors = F)
  genus_kraken <- as.data.frame(table(genusDf$Kraken, useNA = "always"), stringsAsFactors = F)
  mergedTableBC <- merge(genus_blast, genus_centrifuge, by = 1, all = TRUE)
  mergedTableKK <- merge(genus_kaiju, genus_kraken, by = 1, all = TRUE)
  mergedTable <- merge(mergedTableBC, mergedTableKK, by = 1, all = TRUE)
  colnames(mergedTable) <- c("Phylum","Blast","Centrifuge", "Kaiju","Kraken")
  
  # keep only genera with results from all classifiers
  mergedTable <- mergedTable[rowSums(is.na(mergedTable)) < 1 ,]
  
  # start converting dataframe to format for ggplot
  rowNames <- mergedTable[,1]
  rowNamesList <- as.vector(rowNames)
  rowNamesList[is.na(rowNamesList)] <- "Unclassified"
  plotDf <- as.data.frame(mergedTable[,-1])
  rownames(plotDf) <- rowNamesList
  plotDf[is.na(plotDf)] <- 0
  
  noUnclassifiedDf <- plotDf[!(row.names(plotDf)) %in% c("Unclassified", "root"), ]
  blastGGdf <- as.data.frame(noUnclassifiedDf$Blast)
  blastGGdf[,2] <- "Blast"
  blastGGdf[,3] <- rownames(noUnclassifiedDf)
  colnames(blastGGdf) <- c("count", "classifier", "genus")
  
  centGGdf <- as.data.frame(noUnclassifiedDf$Centrifuge)
  centGGdf[,2] <- "Centrifuge"
  centGGdf[,3] <- rownames(noUnclassifiedDf)
  colnames(centGGdf) <- c("count", "classifier", "genus")
  
  kaijuGGdf <- as.data.frame(noUnclassifiedDf$Kaiju)
  kaijuGGdf[,2] <- "Kaiju"
  kaijuGGdf[,3] <- rownames(noUnclassifiedDf)
  colnames(kaijuGGdf) <- c("count", "classifier", "genus")
  
  krakenGGdf <- as.data.frame(noUnclassifiedDf$Kraken)
  krakenGGdf[,2] <- "Kraken2"
  krakenGGdf[,3] <- rownames(noUnclassifiedDf)
  colnames(krakenGGdf) <- c("count", "classifier", "genus")
  
  mergedGGBC <- rbind(blastGGdf, centGGdf)
  mergedGGKK <- rbind(kaijuGGdf, krakenGGdf)
  phylumMergedGGdf <- rbind(mergedGGBC, mergedGGKK)

  phylumMergedGGdf_gt5 <- phylumMergedGGdf[phylumMergedGGdf$count >= 5, ]
  rownames(phylumMergedGGdf_gt5) <- c()
  }

phylumNum <- ggplot(data = phylumMergedGGdf) + 
  geom_bar(data = phylumMergedGGdf, aes(y = count, x = classifier, fill = genus), position="stack", stat="identity", width = 0.5) + 
  scale_fill_manual(values=rep(brewer.pal(12,"Paired"),times=5)) + 
  labs(x="", y="Num. of contigs", fill = "") +
  theme_classic(base_size = 25)
phylumNum

phylumProp <- ggplot(phylumMergedGGdf)+ 
  geom_bar(data = phylumMergedGGdf, aes(y = count, x = classifier, fill = genus), position="fill", stat="identity") + 
  scale_fill_manual(values=rep(brewer.pal(12,"Paired"),times=5)) + 
  labs(x="", y="Prop. annotated contigs", fill = "") +
  theme_classic(base_size = 25)
phylumProp

##### Family
{
  genusDf <- read.table("~/Box/cyanobloom/classifier_family.table", header = TRUE, row.names = 1 , sep = "\t")
  colnames(genusDf) <- c("Blast","Centrifuge", "Kaiju","Kraken","Agreement")
  genus_blast <- as.data.frame(table(genusDf$Blast, useNA = "always"), stringsAsFactors = F)
  genus_centrifuge <- as.data.frame(table(genusDf$Centrifuge, useNA = "always"), stringsAsFactors = F)
  genus_kaiju <- as.data.frame(table(genusDf$Kaiju, useNA = "always"), stringsAsFactors = F)
  genus_kraken <- as.data.frame(table(genusDf$Kraken, useNA = "always"), stringsAsFactors = F)
  mergedTableBC <- merge(genus_blast, genus_centrifuge, by = 1, all = TRUE)
  mergedTableKK <- merge(genus_kaiju, genus_kraken, by = 1, all = TRUE)
  mergedTable <- merge(mergedTableBC, mergedTableKK, by = 1, all = TRUE)
  colnames(mergedTable) <- c("Family","Blast","Centrifuge", "Kaiju","Kraken")
  
  # keep only genera with results from all classifiers
  mergedTable <- mergedTable[rowSums(is.na(mergedTable)) < 1 ,]
  # combine Anabaena and Dolichospermum
  mergedTable[mergedTable$Family == "Aphanizomenonaceae", 2:5] + mergedTable[mergedTable$Family == "Nostocaceae", 2:5]
  # manually insert combined Anabaena/Dolichospermum into Anabaena row
  rownames(mergedTable) <- c()
  mergedTable[mergedTable$Family == "Aphanizomenonaceae", ]
  mergedTable[2,] <- c("Aphanizomenonaceae/Nostocaceae", 775, 766, 608, 762)
  mergedTable <- mergedTable[mergedTable$Family != "Nostocaceae" ,]
  mergedTable <- transform(mergedTable, Blast = as.numeric(Blast), Centrifuge = as.numeric(Centrifuge), Kaiju = as.numeric(Kaiju), Kraken = as.numeric(Kraken))
  #str(mergedTable)
  View(mergedTable)
  
  
  # start converting dataframe to format for ggplot
  rowNames <- mergedTable[,1]
  rowNamesList <- as.vector(rowNames)
  rowNamesList[is.na(rowNamesList)] <- "Unclassified"
  plotDf <- as.data.frame(mergedTable[,-1])
  rownames(plotDf) <- rowNamesList
  plotDf[is.na(plotDf)] <- 0
  
  noUnclassifiedDf <- plotDf[!(row.names(plotDf)) %in% c("Unclassified", "root"), ]
  blastGGdf <- as.data.frame(noUnclassifiedDf$Blast)
  blastGGdf[,2] <- "Blast"
  blastGGdf[,3] <- rownames(noUnclassifiedDf)
  colnames(blastGGdf) <- c("count", "classifier", "genus")
  
  centGGdf <- as.data.frame(noUnclassifiedDf$Centrifuge)
  centGGdf[,2] <- "Centrifuge"
  centGGdf[,3] <- rownames(noUnclassifiedDf)
  colnames(centGGdf) <- c("count", "classifier", "genus")
  
  kaijuGGdf <- as.data.frame(noUnclassifiedDf$Kaiju)
  kaijuGGdf[,2] <- "Kaiju"
  kaijuGGdf[,3] <- rownames(noUnclassifiedDf)
  colnames(kaijuGGdf) <- c("count", "classifier", "genus")
  
  krakenGGdf <- as.data.frame(noUnclassifiedDf$Kraken)
  krakenGGdf[,2] <- "Kraken2"
  krakenGGdf[,3] <- rownames(noUnclassifiedDf)
  colnames(krakenGGdf) <- c("count", "classifier", "genus")
  
  mergedGGBC <- rbind(blastGGdf, centGGdf)
  mergedGGKK <- rbind(kaijuGGdf, krakenGGdf)
  familyMergedGGdf <- rbind(mergedGGBC, mergedGGKK)
  
  familyMergedGGdf_gt5 <- familyMergedGGdf[familyMergedGGdf$count >= 5, ]
  rownames(familyMergedGGdf_gt5) <- c()
}

familyNum <- ggplot(data = familyMergedGGdf_gt5) + 
  geom_bar(data = familyMergedGGdf_gt5, aes(y = count, x = classifier, fill = genus), position="stack", stat="identity", width = 0.5) + 
  scale_fill_manual(values=rep(brewer.pal(12,"Paired"),times=5)) + 
  labs(x="", y="Num. of contigs", fill = "") +
  theme_classic(base_size = 25)
  #theme(axis.text.x=element_text(angle=45,hjust=1))
familyNum

familyProp <- ggplot(familyMergedGGdf_gt5)+ 
  geom_bar(data = familyMergedGGdf_gt5, aes(y = count, x = classifier, fill = genus), position="fill", stat="identity") + 
  scale_fill_manual(values=rep(brewer.pal(12,"Paired"),times=5)) + 
  labs(x="", y="Prop. annotated contigs", fill = "") +
  theme_classic(base_size = 25)
  #theme(axis.text.x=element_text(angle=45,hjust=1))
familyProp

##### Genus 
{
genusDf <- read.table("~/Box/cyanobloom/classifier_genus.table", header = TRUE, row.names = 1 , sep = "\t")
colnames(genusDf) <- c("Blast","Centrifuge", "Kaiju","Kraken","Agreement")
genus_blast <- as.data.frame(table(genusDf$Blast, useNA = "always"), stringsAsFactors = F)
genus_centrifuge <- as.data.frame(table(genusDf$Centrifuge, useNA = "always"), stringsAsFactors = F)
genus_kaiju <- as.data.frame(table(genusDf$Kaiju, useNA = "always"), stringsAsFactors = F)
genus_kraken <- as.data.frame(table(genusDf$Kraken, useNA = "always"), stringsAsFactors = F)
mergedTableBC <- merge(genus_blast, genus_centrifuge, by = 1, all = TRUE)
mergedTableKK <- merge(genus_kaiju, genus_kraken, by = 1, all = TRUE)
mergedTable <- merge(mergedTableBC, mergedTableKK, by = 1, all = TRUE)
colnames(mergedTable) <- c("Genus","Blast","Centrifuge", "Kaiju","Kraken")
#View(mergedTable)

# keep only genera with results from all classifiers
mergedTable <- mergedTable[rowSums(is.na(mergedTable)) < 1 ,]
write(as.vector(mergedTable$Genus), file = "shared_genera.txt")
# combine Anabaena and Dolichospermum
mergedTable[mergedTable$Genus == "Anabaena", 2:5] + mergedTable[mergedTable$Genus == "Dolichospermum", 2:5]

# manually insert combined Anabaena/Dolichospermum into Anabaena row
rownames(mergedTable) <- c()
mergedTable[mergedTable$Genus == "Anabaena", ]
mergedTable[3,] <- c("Anabaena/Dolichospermum", 693, 713, 553, 698)
mergedTable <- mergedTable[mergedTable$Genus != "Dolichospermum" ,]
mergedTable <- transform(mergedTable, Blast = as.numeric(Blast), Centrifuge = as.numeric(Centrifuge), Kaiju = as.numeric(Kaiju), Kraken = as.numeric(Kraken))
mergedTable[mergedTable$Genus == "Nanopelagicus", ]
mergedTable[32,1] <- "Candidatus Nanopelagicus"
mergedTable[mergedTable$Genus == "Planktophila", ]
mergedTable[38,1] <- "Candidatus Planktophila"

#str(mergedTable)
View(mergedTable)
# start converting dataframe to format for ggplot
rowNames <- mergedTable[,1]
rowNamesList <- as.vector(rowNames)
rowNamesList[is.na(rowNamesList)] <- "Unclassified"
plotDf <- as.data.frame(mergedTable[,-1])
rownames(plotDf) <- rowNamesList
plotDf[is.na(plotDf)] <- 0
#View(plotDf)

noUnclassifiedDf <- plotDf[!(row.names(plotDf)) %in% c("Unclassified", "root"), ]
blastGGdf <- as.data.frame(noUnclassifiedDf$Blast)
blastGGdf[,2] <- "Blast"
blastGGdf[,3] <- rownames(noUnclassifiedDf)
colnames(blastGGdf) <- c("count", "classifier", "genus")

centGGdf <- as.data.frame(noUnclassifiedDf$Centrifuge)
centGGdf[,2] <- "Centrifuge"
centGGdf[,3] <- rownames(noUnclassifiedDf)
colnames(centGGdf) <- c("count", "classifier", "genus")

kaijuGGdf <- as.data.frame(noUnclassifiedDf$Kaiju)
kaijuGGdf[,2] <- "Kaiju"
kaijuGGdf[,3] <- rownames(noUnclassifiedDf)
colnames(kaijuGGdf) <- c("count", "classifier", "genus")

krakenGGdf <- as.data.frame(noUnclassifiedDf$Kraken)
krakenGGdf[,2] <- "Kraken2"
krakenGGdf[,3] <- rownames(noUnclassifiedDf)
colnames(krakenGGdf) <- c("count", "classifier", "genus")

mergedGGBC <- rbind(blastGGdf, centGGdf)
mergedGGKK <- rbind(kaijuGGdf, krakenGGdf)
genusMergedGGdf <- rbind(mergedGGBC, mergedGGKK)
#View(mergedGGdf)
rownames(genusMergedGGdf) <- c()
table(genusMergedGGdf$count)

#genusMergedGGdf_gt5 <- genusMergedGGdf[genusMergedGGdf$count >= 5, ]
lessThan5 <- genusMergedGGdf[genusMergedGGdf$count < 5, ]
lessThan5Genera <- table(lessThan5$genus) == 4
lessThan5Genera <- as.data.frame(lessThan5Genera)
lessThan5Genera <- subset(lessThan5Genera, lessThan5Genera == "TRUE" )
lessThan5Genera <- rownames(lessThan5Genera)
lessThan5Genera
}

genusMergedGGdf[genusMergedGGdf$genus %in% lessThan5Genera ] <- "Other"
genusMergedGGdf_gt5 <- genusMergedGGdf
#View(mergedGGdf_gt5)

genusMergedGGdf_gt5$genus <- factor(genusMergedGGdf_gt5$genus, levels = c("Anabaena/Dolichospermum","Aphanizomenon","Calothrix","Cyanobium","Cylindrospermum",
                                                                          "Microcystis","Nostoc","Pseudanabaena","Sphaerospermopsis","Synechococcus",
                                                                          "Acidovorax","Altererythrobacter", "Azospirillum","Bosea","Bradyrhizobium",
                                                                          "Burkholderia","Chelatococcus","Curvibacter","Erythrobacter","Hydrogenophaga",
                                                                          "Limnohabitans","Mesorhizobium","Methylobacterium","Novosphingobium","Phenylobacterium",
                                                                          "Pseudomonas","Rhizobium","Sphingobium","Sphingomonas","Sphingopyxis",
                                                                          "Sphingorhabdus","Variovorax","Brachionus","Candidatus Nanopelagicus",
                                                                          "Candidatus Planktophila","Durinskia","Gossypium", "Other"))

genus_palette <- c(# cyanobacteria
                    "#54de91",
                    "#68743b",
                    "#67e746",
                    "#205f25",
                    "#47b100",
                    "#a1d2aa",
                    "#00c34d",
                    "#9da86e",
                    "#628a00",
                    "#697c27",
                    # proteobacteria
                    "#0045ad",
                     "#7ad8ff",
                     "#ae33b9",
                     "#609dab",
                     "#fb75f8",
                     "#29464d",
                     "#bd68f1",
                     "#00364f",
                     "#ac83ff",
                     "#005c81",
                     "#cd68c4",
                     "#01b1fb",
                     "#500169",
                     "#b2b8c7",
                     "#7f1b7a",
                     "#e5beff",
                     "#25295c",
                     "#be9dff",
                     "#603a5c",
                     "#027df3",
                     "#8c7989",
                     "#bc9bb7",
                   #
                    "yellow",
                   "dark gray",
                   "gray",
                   "red",
                   "orange",
                   "black")
                   
  

genusNum <- ggplot(data = genusMergedGGdf_gt5) + 
  geom_bar(data = genusMergedGGdf_gt5, aes(y = count, x = classifier, fill = genus), position="stack", stat="identity", width = 0.5) + 
  scale_fill_manual(values = genus_palette) + 
  labs(x="", y="Num. of contigs", fill = "") +
  theme_classic(base_size = 25)
genusNum

genusProp <- ggplot(genusMergedGGdf_gt5)+ 
  geom_bar(data = genusMergedGGdf_gt5, aes(y = count, x = classifier, fill = genus), position="fill", stat="identity") + 
  scale_fill_manual(values=genus_palette) + 
  labs(x="", y="Prop. annotated contigs", fill = "") +
  theme_classic(base_size = 25)
genusProp

## Plotting together
gSuperkingdom <- ggplotGrob(superkingdomProp)
gGenus <- ggplotGrob(genusProp)
maxWidth = grid::unit.pmax(gSuperkingdom$widths[2:5], gGenus$widths[2:5])
gSuperkingdom$widths[2:5] <- as.list(maxWidth)
gGenus$widths[2:5] <- as.list(maxWidth)
grid.arrange(gSuperkingdom, gGenus, layout_matrix = rbind(c(1,NA),2), nrow = 2, heights = c(1,2), widths = c(1.5,1))


grid.arrange(superkingdomNum, phylumNum, familyNum, genusNum, nrow = 4)
grid.arrange(phylumProp, familyProp, layout_matrix = rbind(c(1,NA),2), nrow = 2, heights = c(1,1), widths = c(1.5,1))
grid.arrange(superkingdomNum, genusNum, nrow = 2)
grid.arrange(superkingdomProp, genusProp, nrow = 2)

###############################################
### Histogram of unclassified contigs by length ### 

contigLengths <- as.data.frame(read.table("~/Box/cyanobloom/sequences/consensus_contigLengths.tsv", sep = "\t", header = FALSE))
blastContigLengths <- as.data.frame(read.table("~/Box/cyanobloom/sequences/consensus_blastUnclassifiedcontigLengths.tsv", sep = "\t", header = FALSE))
centContigLengths <- as.data.frame(read.table("~/Box/cyanobloom/sequences/consensus_centrifugeUnclassifiedcontigLengths.tsv", sep = "\t", header = FALSE))
kaijuContigLengths <- as.data.frame(read.table("~/Box/cyanobloom/sequences/consensus_kaijuUnclassifiedcontigLengths.tsv", sep = "\t", header = FALSE))
krakenContigLengths <- as.data.frame(read.table("~/Box/cyanobloom/sequences/consensus_krakenUnclassifiedcontigLengths.tsv", sep = "\t", header = FALSE))


hist4 <- ggplot(contigLengths, aes(x = as.numeric(contigLengths$V2))) + 
  geom_histogram(binwidth = 10000, color = "black", fill = "white", show.legend = TRUE) + 
  lims(x = c(0,750000), y = c(0,425)) +
  #geom_histogram(data = blastContigLengths, aes(x= as.numeric(blastContigLengths$V2)), binwidth = 10000, alpha = 1, color = "black", fill = "blue", show.legend = TRUE) +
  #geom_histogram(data = kaijuContigLengths, aes(x= as.numeric(kaijuContigLengths$V2)), binwidth = 10000, alpha = 1, color = "black", fill = "green", show.legend =  TRUE) +
  #geom_histogram(data = krakenContigLengths, aes(x= as.numeric(krakenContigLengths$V2)), binwidth = 10000, alpha = 1, color = "black", fill = "yellow", show.legend =  TRUE) +
  geom_histogram(data = centContigLengths, aes(x= as.numeric(centContigLengths$V2)), binwidth = 10000, alpha = 1, color = "black", fill = "red", show.legend =  TRUE) +
  scale_colour_manual(name="Unclassified Contigs", values=c("blue","red","green", "yellow"),labels=c("Blast","Centrifuge","Kaiju", "Kraken")) +
  labs(x = "Contig Length (bp)", y = "Num. of contigs") +
  theme(text = element_text(size=25))
hist1 + scale_y_continuous(trans='log10')
hist2 + scale_y_continuous(trans='log10')
hist3 + scale_y_continuous(trans='log10')
hist4 + scale_y_continuous(trans='log10')
grid.arrange(hist1 + scale_y_continuous(trans='log10'),
             hist4 + scale_y_continuous(trans='log10'),
             hist2 + scale_y_continuous(trans='log10'),
             hist3 + scale_y_continuous(trans='log10'),
             nrow = 2)
boxplot(contigLengths$V2, blastContigLengths$V2)
boxplot(blastContigLengths$V2)

################################################
### rRNA annotations

rrna <- read.table("~/Desktop/centrifuge_genus.tsv", sep = "\t", header = FALSE)
table_rrna <- as.data.frame(table(rrna$V2))
col3 <- rep("centrifuge", 3)
table_rrna <- cbind(table_rrna, col3)
table_rrna_noroot <- table_rrna[!table_rrna[,1] == "root" , ]

table_rrna_ordered <- table_rrna_noroot[c("Freq","Var1")]
centGGdf_ordered <- centGGdf[c("count","genus")] 
table_rrna_ordered <- cbind(table_rrna_ordered, rep("rRNA", 23))
centGGdf_ordered <- cbind(centGGdf_ordered, rep("DNA", 723))
colnames(centGGdf_ordered) <- c("count", "taxon", "source")
colnames(table_rrna_ordered) <- c("count", "taxon", "source")
centGGdf_gt7 <- centGGdf_ordered[centGGdf_ordered$count >= 7, ]
table_centrifuge <- rbind(table_rrna_ordered, centGGdf_gt7)

rrna_plot <- ggplot(table_centrifuge) +
  geom_bar(aes_(y=table_centrifuge$count, x=table_centrifuge$source, fill=table_centrifuge$taxon), position="fill", stat="identity", width = 0.5) + 
  #scale_fill_manual(values=rep(brewer.pal(12,"Paired"),times=80)) + 
  labs(x="", y="% annotated contigs", fill = "")
rrna_plot + theme_bw()

commonTaxa <- intersect(table_rrna_ordered$taxon, centGGdf_ordered$taxon)
table_centrifuge_commonTaxa <- table_centrifuge[table_centrifuge$taxon %in% commonTaxa , ]

rrna_plot <- ggplot(table_centrifuge_commonTaxa) +
  geom_bar(aes_(y=table_centrifuge_commonTaxa$count, x=table_centrifuge_commonTaxa$source, fill=table_centrifuge_commonTaxa$taxon), position="fill", stat="identity", width = 0.5) + 
  scale_fill_manual(values=rep(brewer.pal(12,"Paired"),times=80)) + 
  ggtitle("Shared Centrifuge Annotations") +
  labs(x="", y="% annotated contigs", fill = "")
rrna_plot + theme_bw()


#########################################################
# Violin Plots comparing data subsets to whole data 
# Superkingdoms
subsets_superkingdom_counts <- read.delim("subsets_superkingdom_blastn_counts.tsv", header = F, sep = "\t", stringsAsFactors = T)
whole_superkingdom_counts <- read.delim("whole_superkingdom_blastn_counts.tsv", header = F, sep = "\t",  stringsAsFactors = T)
superkingdom_count_plot <- ggplot(data = subsets_superkingdom_counts, aes(x=V1, y=V3)) + 
  geom_violin(fill="slateblue", alpha=0.5, trim = F, scale = "width") + 
  xlab("Superkingdom") +
  ylab("% of contigs") +
  geom_point(data = whole_superkingdom_counts, color = "red", size = 3) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(subsets_superkingdom_counts$V1))) +
  theme(text = element_text(size=20))

# Phyla with >1% representation
subsets_phylum_counts <- read.delim("subsets_phylum_blastn_top1perc.tsv", header = F, sep = "\t",  stringsAsFactors = T)
whole_phylum_counts <- read.delim("whole_phylum_blastn_top1perc.tsv", header = F, sep = "\t",  stringsAsFactors = T)
phylum_count_plot <- ggplot(data = subsets_phylum_counts, aes(x=V1, y=V3)) + 
  geom_violin(fill="slateblue", alpha=0.5, trim = F, scale = "width") + 
  xlab("Phylum") +
  ylab("% of contigs") +
  geom_point(data = whole_phylum_counts, color = "red", size = 3) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(subsets_phylum_counts$V1))) +
  theme(text = element_text(size=20)) +
  scale_y_continuous(limits = c(0,90))

# Orders with >1% representation
subsets_order_counts <- read.delim("subsets_order_blastn_top1perc.tsv", header = F, sep = "\t",  stringsAsFactors = T)
whole_order_counts <- read.delim("whole_order_blastn_top1perc.tsv", header = F, sep = "\t",  stringsAsFactors = T)
order_count_plot <- ggplot(data = subsets_order_counts, aes(x=V1, y=V3)) + 
  geom_violin(fill="slateblue", alpha=0.5, trim = F, scale = "width") + 
  xlab("Order") +
  ylab("% of contigs") +
  geom_point(data = whole_order_counts, color = "red", size = 3) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(subsets_order_counts$V1))) +
  theme(text = element_text(size=20)) +
  scale_y_continuous(limits = c(0,90))

# Families with >1% representation
subsets_family_counts <- read.delim("subsets_family_blastn_top1perc.tsv", header = F, sep = "\t",  stringsAsFactors = T)
whole_family_counts <- read.delim("whole_family_blastn_top1perc.tsv", header = F, sep = "\t",  stringsAsFactors = T)
family_count_plot <- ggplot(data = subsets_family_counts, aes(x=V1, y=V3)) + 
  geom_violin(fill="slateblue", alpha=0.5, trim = F, scale = "width") + 
  xlab("Family") +
  ylab("% of contigs") +
  geom_point(data = whole_family_counts, color = "red", size = 3) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(subsets_family_counts$V1))) +
  theme(text = element_text(size=20)) +
  scale_y_continuous(limits = c(0,50))

# Genera with >1% representation
subsets_genus_counts <- read.delim("subsets_genus_blastn_top1perc.tsv", header = F, sep = "\t",  stringsAsFactors = T)
whole_genus_counts <- read.delim("whole_genus_blastn_top1perc.tsv", header = F, sep = "\t",  stringsAsFactors = T)

# Combine Anabaena and Dolichospermum in whole dataset
whole_genus_counts[1,2:3] <- whole_genus_counts[1,2:3] + whole_genus_counts[9,2:3]
whole_genus_counts <- whole_genus_counts[whole_genus_counts$V1 != "Dolichospermum", ]
rownames(whole_genus_counts) <- c()
levels(whole_genus_counts$V1) <- c(levels(whole_genus_counts$V1), "Anabaena/Dolichospermum")
whole_genus_counts[1,1] <- "Anabaena/Dolichospermum"
whole_genus_counts$V1 <- droplevels(whole_genus_counts$V1)
whole_genus_counts$V1 <- factor(whole_genus_counts$V1, levels = c("Anabaena/Dolichospermum","Azospirillum", "Bosea","Candidatus Fonsibacter","Candidatus Nanopelagicus","Candidatus Planktophila", "Cyanobium","Cylindrospermum","Hypericibacter","Inhella","Limnohabitans","Nostoc","Pseudanabaena","Sphaerospermopsis","Sphingomonas","Sphingopyxis","Sphingorhabdus","Stella","Synechococcus","Unannotated"))

# Combine Anabaena and Dolichospermum in subsets
subsets_genus_counts[1,2:3] <- subsets_genus_counts[1,2:3] + subsets_genus_counts[69,2:3]
subsets_genus_counts[2,2:3] <- subsets_genus_counts[2,2:3] + subsets_genus_counts[70,2:3]
subsets_genus_counts[3,2:3] <- subsets_genus_counts[3,2:3] + subsets_genus_counts[71,2:3]
subsets_genus_counts[4,2:3] <- subsets_genus_counts[4,2:3] + subsets_genus_counts[72,2:3]
subsets_genus_counts[5,2:3] <- subsets_genus_counts[5,2:3] + subsets_genus_counts[73,2:3]
subsets_genus_counts[6,2:3] <- subsets_genus_counts[6,2:3] + subsets_genus_counts[74,2:3]
subsets_genus_counts[7,2:3] <- subsets_genus_counts[7,2:3] + subsets_genus_counts[75,2:3]
subsets_genus_counts[8,2:3] <- subsets_genus_counts[8,2:3] + subsets_genus_counts[76,2:3]
subsets_genus_counts[9,2:3] <- subsets_genus_counts[9,2:3] + subsets_genus_counts[77,2:3]
subsets_genus_counts[10,2:3] <- subsets_genus_counts[10,2:3] + subsets_genus_counts[78,2:3]

subsets_genus_counts <- subsets_genus_counts[subsets_genus_counts$V1 != "Dolichospermum", ]
rownames(subsets_genus_counts) <- c()
levels(subsets_genus_counts$V1) <- c(levels(subsets_genus_counts$V1), "Anabaena/Dolichospermum")
subsets_genus_counts[1,1] <- "Anabaena/Dolichospermum" 
subsets_genus_counts[2,1] <- "Anabaena/Dolichospermum" 
subsets_genus_counts[3,1] <- "Anabaena/Dolichospermum" 
subsets_genus_counts[4,1] <- "Anabaena/Dolichospermum" 
subsets_genus_counts[5,1] <- "Anabaena/Dolichospermum" 
subsets_genus_counts[6,1] <- "Anabaena/Dolichospermum" 
subsets_genus_counts[7,1] <- "Anabaena/Dolichospermum" 
subsets_genus_counts[8,1] <- "Anabaena/Dolichospermum" 
subsets_genus_counts[9,1] <- "Anabaena/Dolichospermum" 
subsets_genus_counts[10,1] <- "Anabaena/Dolichospermum" 
subsets_genus_counts$V1 <- droplevels(subsets_genus_counts$V1)
subsets_genus_counts$V1 <- factor(subsets_genus_counts$V1, levels = c("Anabaena/Dolichospermum","Azospirillum", "Bosea","Candidatus Fonsibacter","Candidatus Nanopelagicus","Candidatus Planktophila", "Cyanobium","Cylindrospermum","Hypericibacter","Inhella","Limnohabitans","Nostoc","Pseudanabaena","Sphaerospermopsis","Sphingomonas","Sphingopyxis","Sphingorhabdus","Stella","Synechococcus","Unannotated"))


genus_count_plot <- ggplot(data = subsets_genus_counts, aes(x=V1, y=V3)) + 
  geom_violin(fill="slateblue", alpha=0.5, trim = F, scale = "width") + 
  xlab("Genus") +
  ylab("% of contigs") +
  geom_point(data = whole_genus_counts, color = "red", size = 3) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(subsets_genus_counts$V1))) +
  theme(text = element_text(size=20)) +
  scale_y_continuous(limits = c(0,80))


### Graph based on assembly size rather than % of contigs

# Superkingdoms
subsets_superkingdom_size <- read.delim("~/Box/cyanobloom/subsets_superkingdom_blastn_size.tsv", header = F, sep = "\t",  stringsAsFactors = T)
whole_superkingdom_size <- read.delim("~/Box/cyanobloom/whole_superkingdom_blastn_size.tsv", header = F, sep = "\t",  stringsAsFactors = T)
superkingdom_size_plot <- ggplot(data = subsets_superkingdom_size, aes(x=V1, y=V3)) + 
  geom_violin(fill="slateblue", alpha=0.5, trim = F, scale = "width") + 
  xlab("Superkingdom") +
  ylab("% of assembly size") +
  geom_point(data = whole_superkingdom_size, color = "black", size = 3) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(subsets_superkingdom_size$V1))) +
  theme(text = element_text(size=25))
superkingdom_size_plot
grid.arrange(superkingdom_count_plot, superkingdom_size_plot)

# Phyla with >1% representation
subsets_phylum_size <- read.delim("subsets_phylum_blastn_size_top1perc.tsv", header = F, sep = "\t",  stringsAsFactors = T)
whole_phylum_size <- read.delim("whole_phylum_blastn_size_top1perc.tsv", header = F, sep = "\t",  stringsAsFactors = T)
phylum_size_plot <- ggplot(data = subsets_phylum_size, aes(x=V1, y=V3)) + 
  geom_violin(fill="slateblue", alpha=0.5, trim = F, scale = "width") + 
  xlab("Phylum") +
  ylab("% of assembly size") +
  geom_point(data = whole_phylum_size, color = "red", size = 3) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(subsets_phylum_size$V1))) +
  theme(text = element_text(size=25)) +
  scale_y_continuous(limits = c(0,90))

grid.arrange(phylum_count_plot, phylum_size_plot)

# Orders with >1% representation
subsets_order_size <- read.delim("subsets_order_blastn_size_top1perc.tsv", header = F, sep = "\t",  stringsAsFactors = T)
whole_order_size <- read.delim("whole_order_blastn_size_top1perc.tsv", header = F, sep = "\t",  stringsAsFactors = T)
order_size_plot <- ggplot(data = subsets_order_size, aes(x=V1, y=V3)) + 
  geom_violin(fill="slateblue", alpha=0.5, trim = F, scale = "width") + 
  xlab("Order") +
  ylab("% of assembly size") +
  geom_point(data = whole_order_size, color = "red", size = 3) +
  coord_flip() +
 scale_x_discrete(limits = rev(levels(subsets_order_size$V1))) +
  theme(text = element_text(size=25)) +
  scale_y_continuous(limits = c(0,80))

grid.arrange(order_count_plot, order_size_plot)

# Families with >1% representation
subsets_family_size <- read.delim("subsets_family_blastn_size_top1perc.tsv", header = F, sep = "\t",  stringsAsFactors = T)
whole_family_size <- read.delim("whole_family_blastn_size_top1perc.tsv", header = F, sep = "\t",  stringsAsFactors = T)
family_size_plot <- ggplot(data = subsets_family_size, aes(x=V1, y=V3)) + 
  geom_violin(fill="slateblue", alpha=0.5, trim = F, scale = "width") + 
  xlab("Family") +
  ylab("% of assembly size") +
  geom_point(data = whole_family_size, color = "red", size = 3) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(subsets_family_size$V1))) +
  theme(text = element_text(size=25)) +
  scale_y_continuous(limits = c(0,50))

grid.arrange(family_count_plot, family_size_plot)

# Genera with >1% representation
subsets_genus_size <- read.delim("~/Box/cyanobloom/subsets_genus_blastn_size_top1perc.tsv", header = F, sep = "\t",  stringsAsFactors = T)
whole_genus_size <- read.delim("~/Box/cyanobloom/whole_genus_blastn_size_top1perc.tsv", header = F, sep = "\t",  stringsAsFactors = T)

# Combine Anabaena and Dolichospermum in whole dataset
whole_genus_size[2,2:3] <- whole_genus_size[2,2:3] + whole_genus_size[11,2:3]
whole_genus_size <- whole_genus_size[whole_genus_size$V1 != "Dolichospermum", ]
rownames(whole_genus_size) <- c()
levels(whole_genus_size$V1) <- c(levels(whole_genus_size$V1), "Anabaena/Dolichospermum")
whole_genus_size[2,1] <- "Anabaena/Dolichospermum"
whole_genus_size$V1 <- droplevels(whole_genus_size$V1)
whole_genus_size$V1 <- factor(whole_genus_size$V1, levels = c("Allochromatium","Anabaena/Dolichospermum","Azospirillum", "Bosea","Bradyrhizobium","Candidatus Fonsibacter","Candidatus Nanopelagicus","Candidatus Paracaedimonas","Candidatus Planktophila", "Cyanobium","Erythrobacter","Hypericibacter","Inhella","Limnohabitans","Magnetospirillum","Pseudanabaena","Raphidiopsis","Rhizobium","Rhodopseudomonas","Rhodospirillum","Sphingobium","Sphingomonas","Sphingopyxis","Sphingorhabdus","Stella","Synechococcus","Unannotated"))

# Combine Anabaena and Dolichospermum in subsets
subsets_genus_size[2,2:3] <- subsets_genus_size[2,2:3] + subsets_genus_size[70,2:3]
subsets_genus_size[3,2:3] <- subsets_genus_size[3,2:3] + subsets_genus_size[71,2:3]
subsets_genus_size[4,2:3] <- subsets_genus_size[4,2:3] + subsets_genus_size[72,2:3]
subsets_genus_size[5,2:3] <- subsets_genus_size[5,2:3] + subsets_genus_size[73,2:3]
subsets_genus_size[6,2:3] <- subsets_genus_size[6,2:3] + subsets_genus_size[74,2:3]
subsets_genus_size[7,2:3] <- subsets_genus_size[7,2:3] + subsets_genus_size[75,2:3]
subsets_genus_size[8,2:3] <- subsets_genus_size[8,2:3] + subsets_genus_size[76,2:3]
subsets_genus_size[9,2:3] <- subsets_genus_size[9,2:3] + subsets_genus_size[77,2:3]
subsets_genus_size[10,2:3] <- subsets_genus_size[10,2:3] + subsets_genus_size[78,2:3]
subsets_genus_size[11,2:3] <- subsets_genus_size[11,2:3] + subsets_genus_size[79,2:3]

subsets_genus_size <- subsets_genus_size[subsets_genus_size$V1 != "Dolichospermum", ]
rownames(subsets_genus_size) <- c()
levels(subsets_genus_size$V1) <- c(levels(subsets_genus_size$V1), "Anabaena/Dolichospermum")
subsets_genus_size[2,1] <- "Anabaena/Dolichospermum" 
subsets_genus_size[3,1] <- "Anabaena/Dolichospermum" 
subsets_genus_size[4,1] <- "Anabaena/Dolichospermum" 
subsets_genus_size[5,1] <- "Anabaena/Dolichospermum" 
subsets_genus_size[6,1] <- "Anabaena/Dolichospermum" 
subsets_genus_size[7,1] <- "Anabaena/Dolichospermum" 
subsets_genus_size[8,1] <- "Anabaena/Dolichospermum" 
subsets_genus_size[9,1] <- "Anabaena/Dolichospermum" 
subsets_genus_size[10,1] <- "Anabaena/Dolichospermum" 
subsets_genus_size[11,1] <- "Anabaena/Dolichospermum" 
subsets_genus_size$V1 <- droplevels(subsets_genus_size$V1)
subsets_genus_size$V1 <- factor(subsets_genus_size$V1, levels = c("Allochromatium","Anabaena/Dolichospermum","Azospirillum", "Bosea","Bradyrhizobium","Candidatus Fonsibacter","Candidatus Nanopelagicus","Candidatus Paracaedimonas","Candidatus Planktophila", "Cyanobium","Erythrobacter","Hypericibacter","Inhella","Limnohabitans","Magnetospirillum","Pseudanabaena","Raphidiopsis","Rhizobium","Rhodopseudomonas","Rhodospirillum","Sphingobium","Sphingomonas","Sphingopyxis","Sphingorhabdus","Stella","Synechococcus","Unannotated"))


genus_size_plot <- ggplot(data = subsets_genus_size, aes(x=V1, y=V3)) + 
  geom_violin(fill="slateblue", alpha=0.5, trim = F, scale = "width") + 
  xlab("Genus") +
  ylab("% of assembly size") +
  geom_point(data = whole_genus_size, color = "black", size = 3) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(subsets_genus_size$V1))) +
  theme(text = element_text(size=25)) +
  scale_y_continuous(limits = c(0,55))
genus_size_plot
grid.arrange(genus_count_plot, genus_size_plot)




#### Plot together
gSuperkingdom <- ggplotGrob(superkingdom_size_plot)
gGenus <- ggplotGrob(genus_size_plot)
maxWidth = grid::unit.pmax(gSuperkingdom$widths[2:5], gGenus$widths[2:5])
gSuperkingdom$widths[2:5] <- as.list(maxWidth)
gGenus$widths[2:5] <- as.list(maxWidth)
grid.arrange(gSuperkingdom, gGenus, nrow = 2, heights = c(1,3))

grid.arrange(superkingdom_size_plot, family_size_plot, phylum_size_plot, genus_size_plot, nrow = 2, widths = c(1,2))

######################################################
# COG summary  of predicted proteins
library(ggplot2)
library(grid)
setwd("~/Box/cyanobloom/sequences")
cog_table <- as.data.frame(read.table("cog_table.tsv", header = T, sep = "\t"))

coglegend <- "A: RNA processing and modification
B: Chromatin structure and dynamics
C: Energy production and conversion
D: Cell cycle control, cell division, chromosome partitioning 
E: Amino acid transport and metabolism
F: Nucleotide transport and metabolism
G: Carbohydrate transport and metabolism
H: Coenzyme transport and metabolism
I: Lipid transport and metabolism
J: Translation, ribosomal structure and biogenesis
K: Transcription
L: Replication, recombination and repair
M: Cell wall/membrane/envelope biogenesis
N: Cell Motility
O: Posttranslational modification, protein turnover, chaperones
P: Inorganic ion transport and metabolism
Q: Secondary metabolites biosynthesis, transport and catabolism
S: Function unknown
T: Signal transduction mechanisms
U: Intracellular trafficking, secretion, and vesicular transport
V: Defense mechanisms
W: Extracellular structures
Y: Nuclear structure
Z: Cytoskeleton"
cog_table
cog_hist <- ggplot(data = cog_table) +
  theme(plot.margin=unit(c(0.5,5,0.5,0.5),"in"),
        panel.background = element_rect("white", "white"),
        panel.grid.major = element_line("gray", 0.25),
        panel.grid.minor = element_line("gray", 0.1)) +
  geom_histogram(aes(COG), stat = "count", show.legend = F) +
  ylab("Count") +
  coord_cartesian(clip = 'off') +
  annotation_custom(grob = textGrob(coglegend, hjust = 0), ymin = 0, ymax = 9000, xmin = 20, xmax = 30)

cog_hist                    
                    
                    
cog_presence <- ggplot(cog_table, aes(COG, Taxon)) + 
  geom_tile() +
  ylab("Genus") +
  theme(plot.margin=unit(c(0.5,5,0.5,0.5),"in"),
        panel.background = element_rect("white", "white"),
        panel.grid.major.x = element_line("gray", 0.2, "dashed"),
        panel.grid.major.y = element_line("gray", 0.2)
        #panel.grid.minor = element_line("gray", 0.1)
        ) +
  coord_cartesian(clip = 'off') +
  annotation_custom(grob = textGrob(coglegend, hjust = 0), ymin = 100, ymax = 100, xmin = 20, xmax = 30) +
  scale_x_discrete(limits = rev(levels(cog_table$Taxon)), )
cog_presence

genus_hist <- ggplot(data = cog_table) +
  theme(plot.margin=unit(c(0.5,5,0.5,0.5),"in"),
        panel.background = element_rect("white", "white"),
        panel.grid.major = element_line("gray", 0.25),
        panel.grid.minor = element_line("gray", 0.1)) +
  geom_histogram(aes(Taxon), stat = "count", show.legend = F) +
  ylab("Count") +
  coord_cartesian(clip = 'off')
  #annotation_custom(grob = textGrob(coglegend, hjust = 0), ymin = 0, ymax = 9000, xmin = 20, xmax = 30) 
genus_hist

library(pheatmap)
library(grid) 
cog_matrix <- read.delim("cog_matrix.tsv", sep = "\t", header = T, row.names = 1)
cog_matrix <- cog_matrix[ order(row.names(cog_matrix)), ]

cog_matrix["Unknown",] <- cog_matrix["root",] + cog_matrix["Unknown",]
cog_matrix <- cog_matrix[rownames(cog_matrix) != "root",]
cog_matrix <- cog_matrix[rowSums(cog_matrix) > 1, ] # remove singletons which is largest number of genera
cog_matrix_high <- cog_matrix[rowSums(cog_matrix) > 100, ] # get only most represented genera

# combine Anabaena and Dolichospermum
cog_matrix_high[rownames(cog_matrix_high) == "Anabaena", ]
cog_matrix_high[rownames(cog_matrix_high) == "Dolichospermum", ]
cog_matrix_high[rownames(cog_matrix_high) == "Anabaena", ] + cog_matrix_high[rownames(cog_matrix_high) == "Dolichospermum", ]

# manually insert combined Anabaena/Dolichospermum into Anabaena row
cog_matrix_high["Anabaena",] <- c(4,1,549,148,455,149,317,344,191,335,431,894,542,93,392,402,297,2170,450,190,199,0,1,4)
cog_matrix_high <- cog_matrix_high[rownames(cog_matrix_high) != "Dolichospermum" ,]
rownames(cog_matrix_high)[rownames(cog_matrix_high) == "Anabaena"] <- "Anabaena/Dolichospermum"
rownames(cog_matrix_high)[rownames(cog_matrix_high) == "Nanopelagicus"] <- "Candidatus Nanopelagicus"
rownames(cog_matrix_high)[rownames(cog_matrix_high) == "Planktophila"] <- "Candidatus Planktophila"

cog_matrix_high <- cog_matrix_high[order(rownames(cog_matrix_high)),]

cog_matrix_high_log <- log(cog_matrix_high + 1)

pheatmap(as.matrix(cog_matrix), cluster_rows = FALSE, cluster_cols = FALSE)
cog_matrix_heatmap <- pheatmap(as.matrix(cog_matrix_high_log), cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 25)

# function to rotate labels in pheatmap
draw_colnames_0 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 1, hjust = 1, rot = 0, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_0",
                  ns=asNamespace("pheatmap"))

# histograms
melt_cog_matrix_high <- melt(cog_matrix_high)
View(melt_cog_matrix_high)
cog_matrix_high_hist <- ggplot(data = melt_cog_matrix_high) +
  geom_bar(data = melt_cog_matrix_high, aes(x = variable, y = value), stat = "identity") +
  ylab("") +
  theme_bw() +
  theme(#plot.margin=unit(c(0.5,5,0.5,0.5),"in"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_rect("white", "white"),
        panel.grid.major = element_line("dark gray", 0.25),
        panel.grid.minor = element_line("gray", 0.1),
        panel.border = element_blank(),
        text = element_text(size=25))
cog_matrix_high_hist

melt_genus_matrix_high <- melt(t(cog_matrix_high))
genus_matrix_high_hist <- ggplot(data = melt_genus_matrix_high) +
  geom_bar(data = melt_genus_matrix_high, aes(x = Var2, y = value), stat = "identity") +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background = element_rect("white", "white"),
        panel.grid.major = element_line("dark gray", 0.25),
        panel.grid.minor = element_line("gray", 0.1),
        panel.border = element_blank(),
        text = element_text(size=25)) +
        coord_flip() +
        scale_y_reverse() +
        scale_x_discrete(limits = rev(levels(melt_genus_matrix_high$Var2)))
genus_matrix_high_hist

grid.arrange(cog_matrix_high_hist, layout_matrix = rbind(c(NA,1),c(NA,NA)), heights = c(1,4), widths = c(1,4))
grid.arrange(genus_matrix_high_hist, cog_matrix_heatmap[[4]], widths = c(1,4))
grid.arrange(grobs = list(cog_matrix_high_hist, genus_matrix_high_hist, cog_matrix_heatmap[[4]]), layout_matrix = rbind(c(NA,1),c(2,3)), widths = c(1,))


########################################
# Genus assembly size

genusAssemblyStats <- read.delim("~/Box/cyanobloom/genus_assembly-stats_all.tsv")
View(genusAssemblyStats)

assemblySizePlot <- ggplot(data=genusAssemblyStats) +
  geom_col(data = genusAssemblyStats, aes(filename, total_length))

assemblySizePlot
plot(genusAssemblyStats$number~genusAssemblyStats$total_length, xlim = c(3.5e+06, 2.0e+07))
text(genusAssemblyStats$number~genusAssemblyStats$total_length, labels= genusAssemblyStats$filename)
summary(lm(genusAssemblyStats$number~genusAssemblyStats$total_length))
cor(genusAssemblyStats$number, genusAssemblyStats$total_length, method = "kendall")
cor.test(genusAssemblyStats$number, genusAssemblyStats$total_length, method="kendall") 


##########################################
# BLASTn alignment lengths
alignlen <- read.delim("~/Box/cyanobloom/blastn.alignmentlenths", header = F)

summary(alignlen$V1)

hitlen <- read.delim("~/Box/cyanobloom/centrifuge.hitlengths", header = F)
summary(as.numeric(hitlen$V1))



###############################################
# Species Richness vs. Amt of Data

speciesRichnessDf <- read.table("~/Box/cyanobloom/species_richness.tsv", header = T, sep = "\t")
blastnSpeciesRichness <- speciesRichnessDf[speciesRichnessDf$Tool == "Blastn", ]
centrifugeSpeciesRichness <- speciesRichnessDf[speciesRichnessDf$Tool == "Centrifuge", ]
kaijuSpeciesRichness <- speciesRichnessDf[speciesRichnessDf$Tool == "Kaiju", ]
krakenSpeciesRichness <- speciesRichnessDf[speciesRichnessDf$Tool == "Kraken", ]

blastAsymp <- blastnSpeciesRichness[blastnSpeciesRichness$GB != "11", ]
centrifugeAsymp <- centrifugeSpeciesRichness[centrifugeSpeciesRichness$GB != "11", ]
kaijuAsymp <- kaijuSpeciesRichness[kaijuSpeciesRichness$GB != "11", ]
krakenAsymp <- krakenSpeciesRichness[krakenSpeciesRichness$GB != "11", ]

blastGenusNLS <- nls(Shared_Genera ~ SSlogis(GB, Asym, xmid, scal), data = blastnSpeciesRichness[blastnSpeciesRichness$GB != "11", ])
centrifugeGenusNLS <- nls(Shared_Genera ~ SSlogis(GB, Asym, xmid, scal), data = centrifugeSpeciesRichness[centrifugeSpeciesRichness$GB != "11", ])
kaijuGenusNLS <- nls(Shared_Genera ~ SSlogis(GB, Asym, xmid, scal), data = kaijuSpeciesRichness[kaijuSpeciesRichness$GB != "11", ])
krakenGenusNLS <- nls(Shared_Genera ~ SSlogis(GB, Asym, xmid, scal), data = krakenSpeciesRichness[krakenSpeciesRichness$GB != "11", ])
blastFamilyNLS <- nls(Shared_Families ~ SSlogis(GB, Asym, xmid, scal), data = blastnSpeciesRichness[blastnSpeciesRichness$GB != "11", ])
centrifugeFamilyNLS <- nls(Shared_Families ~ SSlogis(GB, Asym, xmid, scal), data = centrifugeSpeciesRichness[centrifugeSpeciesRichness$GB != "11", ])
kaijuFamilyNLS <- nls(Shared_Families ~ SSlogis(GB, Asym, xmid, scal), data = kaijuSpeciesRichness[kaijuSpeciesRichness$GB != "11", ])
krakenFamilyNLS <- nls(Shared_Families ~ SSlogis(GB, Asym, xmid, scal), data = krakenSpeciesRichness[krakenSpeciesRichness$GB != "11", ])
blastOrderNLS <- nls(Shared_Orders ~ SSlogis(GB, Asym, xmid, scal), data = blastnSpeciesRichness[blastnSpeciesRichness$GB != "11", ])
centrifugeOrderNLS <- nls(Shared_Orders ~ SSlogis(GB, Asym, xmid, scal), data = centrifugeSpeciesRichness[centrifugeSpeciesRichness$GB != "11", ])
kaijuOrderNLS <- nls(Shared_Orders ~ SSlogis(GB, Asym, xmid, scal), data = kaijuSpeciesRichness[kaijuSpeciesRichness$GB != "11", ])
krakenOrderNLS <- nls(Shared_Orders ~ SSlogis(GB, Asym, xmid, scal), data = krakenSpeciesRichness[krakenSpeciesRichness$GB != "11", ])
blastPhylumNLS <- nls(Shared_Phyla ~ SSlogis(GB, Asym, xmid, scal), data = blastnSpeciesRichness[blastnSpeciesRichness$GB != "11", ])
centrifugePhylumNLS <- nls(Shared_Phyla ~ SSlogis(GB, Asym, xmid, scal), data = centrifugeSpeciesRichness[centrifugeSpeciesRichness$GB != "11", ])
kaijuPhylumNLS <- nls(Shared_Phyla ~ SSlogis(GB, Asym, xmid, scal), data = kaijuSpeciesRichness[kaijuSpeciesRichness$GB != "11", ])
krakenPhylumNLS <- nls(Shared_Phyla ~ SSlogis(GB, Asym, xmid, scal), data = krakenSpeciesRichness[krakenSpeciesRichness$GB != "11", ])

blastGenusLog = lm(blastnSpeciesRichness$Shared_Genera ~ log(blastnSpeciesRichness$GB))
plot(blastnSpeciesRichness$Shared_Genera ~ blastnSpeciesRichness$GB)
lines(blastnSpeciesRichness$GB, predict(blastGenusLog))


#### Plotting
par(mfrow = c(2,2))
par(xpd=TRUE)
# Genus
{
plot(speciesRichnessDf$Shared_Genera ~ speciesRichnessDf$GB, xlab = "Data Amount (Gbp)", ylab = "Genus Richness", pch = 20)
points(blastnSpeciesRichness$Shared_Genera ~ blastnSpeciesRichness$GB, col = "black", pch = 15)
points(centrifugeSpeciesRichness$Shared_Genera ~ centrifugeSpeciesRichness$GB, col = "blue", pch = 16)
points(kaijuSpeciesRichness$Shared_Genera ~ kaijuSpeciesRichness$GB, col = "red", pch = 17)
points(krakenSpeciesRichness$Shared_Genera ~ krakenSpeciesRichness$GB, col = "green", pch = 18)
lines(blastAsymp$GB, predict(blastGenusNLS), color = "black", lwd = 2)
lines(centrifugeAsymp$GB, predict(centrifugeGenusNLS), col = "blue", lwd = 2)
lines(kaijuAsymp$GB, predict(kaijuGenusNLS), col = "red", lwd = 2)
lines(krakenAsymp$GB, predict(krakenGenusNLS), col = "green", lwd = 2)
}
# Family
{
plot(speciesRichnessDf$Shared_Families ~ speciesRichnessDf$GB, xlab = "Data Amount (Gbp)", ylab = "Family Richness", pch = 20)
points(blastnSpeciesRichness$Shared_Families ~ blastnSpeciesRichness$GB, col = "black", pch = 15)
points(centrifugeSpeciesRichness$Shared_Families ~ centrifugeSpeciesRichness$GB, col = "blue", pch = 16)
points(kaijuSpeciesRichness$Shared_Families ~ kaijuSpeciesRichness$GB, col = "red", pch = 17)
points(krakenSpeciesRichness$Shared_Families ~ krakenSpeciesRichness$GB, col = "green", pch = 18)
lines(blastAsymp$GB, predict(blastFamilyNLS), color = "black", lwd = 2)
lines(centrifugeAsymp$GB, predict(centrifugeFamilyNLS), col = "blue", lwd = 2)
lines(kaijuAsymp$GB, predict(kaijuFamilyNLS), col = "red", lwd = 2)
lines(krakenAsymp$GB, predict(krakenFamilyNLS), col = "green", lwd = 2)
}
# Order
{
plot(speciesRichnessDf$Shared_Orders ~ speciesRichnessDf$GB, xlab = "Data Amount (Gbp)", ylab = "Order Richness", pch = 20)
points(blastnSpeciesRichness$Shared_Orders ~ blastnSpeciesRichness$GB, col = "black", pch = 15)
points(centrifugeSpeciesRichness$Shared_Orders ~ centrifugeSpeciesRichness$GB, col = "blue", pch = 16)
points(kaijuSpeciesRichness$Shared_Orders ~ kaijuSpeciesRichness$GB, col = "red", pch = 17)
points(krakenSpeciesRichness$Shared_Orders ~ krakenSpeciesRichness$GB, col = "green", pch = 18)
lines(blastAsymp$GB, predict(blastOrderNLS), color = "black", lwd = 2)
lines(centrifugeAsymp$GB, predict(centrifugeOrderNLS), col = "blue", lwd = 2)
lines(kaijuAsymp$GB, predict(kaijuOrderNLS), col = "red", lwd = 2)
lines(krakenAsymp$GB, predict(krakenOrderNLS), col = "green", lwd = 2)
}
# Phylum
{
plot(speciesRichnessDf$Shared_Phyla ~ speciesRichnessDf$GB, xlab = "Data Amount (Gbp)", ylab = "Phylum Richness", pch = 20)
points(blastnSpeciesRichness$Shared_Phyla ~ blastnSpeciesRichness$GB, col = "black", pch = 15)
points(centrifugeSpeciesRichness$Shared_Phyla ~ centrifugeSpeciesRichness$GB, col = "blue", pch = 16)
points(kaijuSpeciesRichness$Shared_Phyla ~ kaijuSpeciesRichness$GB, col = "red", pch = 17)
points(krakenSpeciesRichness$Shared_Phyla ~ krakenSpeciesRichness$GB, col = "green", pch = 18)
lines(blastAsymp$GB, predict(blastPhylumNLS), color = "black", lwd = 2)
lines(centrifugeAsymp$GB, predict(centrifugePhylumNLS), col = "blue", lwd = 2)
lines(kaijuAsymp$GB, predict(kaijuPhylumNLS), col = "red", lwd = 2)
lines(krakenAsymp$GB, predict(krakenPhylumNLS), col = "green", lwd = 2)
legend(1.5,30, 
       legend = c("Blastn", "Centrifuge", "Kaiju", "Kraken"), 
       col = c("black","blue","red","green"), 
       pch = c(15,16,17,18), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = T)
}


###### Metagenome Contig Coverage
contigCoverage<- read.table("~/Box/cyanobloom/sequences/assembly_info_coverage.txt", header = T, sep = "\t")
cc <- ggplot(data = contigCoverage) +
      geom_histogram(data = contigCoverage, )

