setwd("~/Desktop")
require(gridExtra)
require(ggplot2)
require(RColorBrewer)
pal <- brewer.pal(11, "Paired")
assemblyDF <- as.data.frame(read.csv("Hornwort_assembly_stats.csv", header = TRUE))
assemblyDF$Sample <- factor(assemblyDF$Sample, levels = assemblyDF$Sample)

plot1leg <- ggplot(data=assemblyDF, aes(x=Sample, y=Genome.Size..Mb., fill = Sample)) + 
  geom_bar(position="dodge", stat="identity") +
  xlab("") +
  ylab("Genome Size (Mb)") +
  scale_fill_manual(values=pal)
plot1leg <- plot1leg + theme_bw() + theme(axis.text.x = element_blank())

### Make separate legend
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend <- get_legend(plot1leg)

plot1 <- ggplot(data=assemblyDF, aes(x=Sample, y=Genome.Size..Mb.)) + 
  geom_bar(position="dodge", stat="identity", fill = pal) +
  xlab("") +
  ylab("Genome Size (Mb)") +
  ggtitle("Genome Size")
plot1 <- plot1 + theme_bw() + theme(axis.text.x = element_blank())

plot2 <- ggplot(data=assemblyDF, aes(x=Sample, y=N50..Mb.)) + 
  geom_bar(position="dodge", stat="identity", fill = pal) +
  xlab("") +
  ylab("N50 (Mb)") +
  ggtitle("N50")
plot2 <- plot2 + theme_bw() + theme(axis.text.x = element_blank())

plot4 <- ggplot(data=assemblyDF, aes(x=Sample, y=N50n)) + 
  geom_bar(position="dodge", stat="identity", fill = pal) +
  xlab("") +
  ylab("N50n") +
  ggtitle("N50n")
plot4 <- plot4 + theme_bw() + theme(axis.text.x = element_blank())

plot3 <- ggplot(data=assemblyDF, aes(x=Sample, y=N90..Mb.)) + 
  geom_bar(position="dodge", stat="identity", fill = pal) +
  xlab("") +
  ylab("N90 (Mb)") +
  ggtitle("N90")
plot3 <- plot3 + theme_bw() + theme(axis.text.x = element_blank())

plot5 <- ggplot(data=assemblyDF, aes(x=Sample, y=N90n)) + 
  geom_bar(position="dodge", stat="identity", fill = pal) +
  xlab("") +
  ylab("N90n") +
  ggtitle("N90n")
plot5 <- plot5 + theme_bw() + theme(axis.text.x = element_blank()) + theme(legend.text=element_text(size=14), legend.title=element_text(size=14))

grid.arrange(plot1, plot2, legend, plot3, layout_matrix = rbind(c(1, 2), c(3, 4)))


##### Proteome

proteomeDF <- as.data.frame(read.csv("Hornwort_proteome_stats.csv", header = TRUE))
proteomeDF$Sample <- factor(proteomeDF$Sample, levels = proteomeDF$Sample)

plot1leg <- ggplot(data=proteomeDF, aes(x=Sample, y=Proteome.Size, fill = Sample)) + 
  geom_bar(position="dodge", stat="identity") +
  xlab("") +
  ylab("Genome Size (Mb)") +
  scale_fill_manual(values=pal)
plot1leg <- plot1leg + theme_bw() + theme(axis.text.x = element_blank())

### Make separate legend
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend <- get_legend(plot1leg)

plot1 <- ggplot(data=proteomeDF, aes(x=Sample, y=Proteome.Size)) + 
  geom_bar(position="dodge", stat="identity", fill = pal) +
  xlab("") +
  ylab("Proteome Size (Ma)") +
  ggtitle("Proteome Size")
plot1 <- plot1 + theme_bw() + theme(axis.text.x = element_blank())

plot2 <- ggplot(data=proteomeDF, aes(x=Sample, y=X..genes)) + 
  geom_bar(position="dodge", stat="identity", fill = pal) +
  xlab("") +
  ylab("No. of genes") +
  ggtitle("No. of genes")
plot2 <- plot2 + theme_bw() + theme(axis.text.x = element_blank())

plot4 <- ggplot(data=proteomeDF, aes(x=Sample, y=Avg..gene.length)) + 
  geom_bar(position="dodge", stat="identity", fill = pal) +
  xlab("") +
  ylab("Avg. Gene Length (bp)") +
  ggtitle("Avg. Gene Length")
plot4 <- plot4 + theme_bw() + theme(axis.text.x = element_blank())

plot3 <- ggplot(data=proteomeDF, aes(x=Sample, y=Avg..exon..length)) + 
  geom_bar(position="dodge", stat="identity", fill = pal) +
  xlab("") +
  ylab("Avg. Exon Length (bp)") +
  ggtitle("Avg. Exon Length")
plot3 <- plot3 + theme_bw() + theme(axis.text.x = element_blank())

grid.arrange(plot2,plot4, legend, plot3, layout_matrix = rbind(c(1, 2), c(3, 4)))


##### PPR
pprDF <- as.data.frame(read.csv("Hornwort_PPR_stats.csv", header = TRUE))
pprDF$Sample <- factor(pprDF$Sample, levels = pprDF$Sample)
pprplot <- ggplot(data=pprDF, aes(x=Sample, y=PPRs, fill = Sample)) + 
  geom_bar(position="dodge", stat="identity") +
  xlab("") +
  ylab("No. PPRs") +
  ggtitle("No. PPRs") +
  scale_fill_manual(values=pal)
pprplot + theme_bw() + theme(axis.text.x = element_blank())



### Gene Presence/Absence Trees
require(ape)
require(phytools)
genomesTree <- read.tree("Hornwort_genomes_ASTRAL.newick")
stomataGenes <- as.matrix(read.csv("CCM_genes.csv", header = TRUE, row.names = 1))

for (i in 1:11){
print(sort(genomesTree$tip.label)[i] == sort(rownames(stomataGenes))[i])
}
  
phylo.heatmap(genomesTree, stomataGenes, fsize = 1, legend = FALSE, mar = c(1,1,5,1), colors = c("gray","white"))
plotTree(genomesTree)











