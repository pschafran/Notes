setwd("/Applications/Phylogenetics/PURC/2017/Combined")
library(ggplot2)
purcdata <- read.csv("ClusterStats_allTotalClusters_2017.csv", head=T)
pcrreps <- read.csv("ClusterStats_PCRreplicates.csv", head=T)
df = data.frame(x = c(purcdata$X990,purcdata$X995,purcdata$default,purcdata$X997,purcdata$expected), type=rep(c("990 all", "995 all", "default", "997 all", "expected")), c(length(purcdata$X990.all),length(purcdata$X995.all),length(purcdata$default),length(purcdata$X997.all),length(purcdata$expected)))

diff = data.frame(x = c(purcdata$expected990Diff,purcdata$expected995Diff,purcdata$expectedDefaultDiff,purcdata$expected997Diff), type=rep(c("990 all", "995 all", "default", "997 all")), c(length(purcdata$X990.all),length(purcdata$X995.all),length(purcdata$default),length(purcdata$X997.all)))

reps = data.frame(x = c(pcrreps$X990,pcrreps$X995,pcrreps$X997,pcrreps$default), type=rep(c("990 all", "995 all", "997 all","Default")), c(length(pcrreps$X990),length(pcrreps$X995),length(pcrreps$X997),length(pcrreps$default)))

ggplot(df, aes(x=df$ClusterStats_purc_run_LFY_990_990_990_990_4_LFY_clustered_reconsensus.fa.2.csv, fill = df$ClusterStats_purc_run_LFY_995_995_995_995_4_LFY_clustered_reconsensus.fa.2.csv)) + geom_histogram(binwidth=1, alpha = 0.5, position = "identity", color = "black", fill = "white")

qbarplot(table(df$ClusterStats_purc_run_LFY_995_995_995_995_4_LFY_clustered_reconsensus.fa.2.csv))

ggplot(df, aes(x, fill = type)) + geom_bar(position = 'dodge') + scale_x_continuous(breaks = seq(1,4,1), limits = c(0,5)) + ylab("Count") + xlab("Alleles per Individual")
ggplot(diff, aes(x, fill = type)) + geom_bar(position = 'dodge') + scale_x_continuous(breaks = seq(-4,4,1), limits = c(-4,4)) + ylab("Count") + xlab("Difference between Expected and Inferred Alleles")
ggplot(reps, aes(x, fill = type)) + geom_bar(position = 'dodge') + scale_x_continuous(breaks = seq(0,4,1), limits = c(-1,5)) + ylab("Count") + xlab("Difference in Inferred Alleles")



###24 April 2018

setwd("/Applications/Phylogenetics/PURC/2018/")
library(ggplot2)
purcdata <- read.csv("ClusterStats_allTotalClusters.csv", head=T)
pcrreps <- read.csv("ClusterStats_replicateSTDEV.csv", head=T)
df = data.frame(x = c(purcdata$X990,purcdata$X995,purcdata$default,purcdata$X997,purcdata$expected), type=rep(c("990 all", "995 all", "default", "997 all", "expected")), c(length(purcdata$X990.all),length(purcdata$X995.all),length(purcdata$default),length(purcdata$X997.all),length(purcdata$expected)))

diff = data.frame(x = c(purcdata$expected990Diff,purcdata$expected995Diff,purcdata$expectedDefaultDiff,purcdata$expected997Diff), type=rep(c("990 all", "995 all", "default", "997 all")), c(length(purcdata$X990.all),length(purcdata$X995.all),length(purcdata$default),length(purcdata$X997.all)))

reps = data.frame(x = c(pcrreps$stdev990,pcrreps$stdev995,pcrreps$stdev997,pcrreps$stdevDefault), type=rep(c("990 all", "995 all", "997 all","Default")), c(length(pcrreps$stdev990),length(pcrreps$stdev995),length(pcrreps$stdev997),length(pcrreps$stdevDefault)))

ggplot(df, aes(x, fill = type)) + geom_bar(position = 'dodge') + scale_x_continuous(breaks = seq(1,4,1), limits = c(0,5)) + ylab("Count") + xlab("Alleles per Individual")
ggplot(diff, aes(x, fill = type)) + geom_bar(position = 'dodge') + scale_x_continuous(breaks = seq(-4,4,1), limits = c(-4,4)) + ylab("Count") + xlab("Difference between Expected and Inferred Alleles")
ggplot(reps, aes(x, fill = type)) + geom_bar(position = 'stack') + scale_x_continuous(breaks = seq(0,4,1), limits = c(-1,5)) + ylab("Count") + xlab("Standard Deviation of Alleles Inferred for Replicates")

###25 Feb 2019
setwd("/Applications/Phylogenetics/PURC/2019/purc_out_BCPrimers_All_10lqbp_997_995_990_997")
flowcyt <- read.csv("FlowCyt.csv")
diploids <- subset(flowcyt, flowcyt$X1 == 1)
tetraploids <- subset(flowcyt, flowcyt$X1 == 2)
hexaploids <- subset(flowcyt, flowcyt$X1 == 3)
octoploids <- subset(flowcyt, flowcyt$X1 == 4)
boxplot(diploids$X2.25, tetraploids$X2.25, hexaploids$X2.25, xlab = "Number of genomes", ylab = "C-values", names = c(1,2,3))

###3 March 2019
###Examining BLASTed LFY hit scores and lengths to determine cutoff for bad hits
setwd("/Applications/Phylogenetics/PURC/2019/ManualDemultiplex")
bitscore <- read.csv("LFY_bitscores.csv", header = FALSE)
summary(bitscore$V1)
boxplot(bitscore$V1)
hist(bitscore$V1)

hitLength <- read.csv("LFY_hitLength.csv", header = FALSE)
summary(hitLength$V1)
boxplot(hitLength$V1)
hist(hitLength$V1)
### Lower cutoff set at 1250 for bit scores and 800 for hit lengths to try to capture outgroups
PGI_BLASTOUT <- read.delim("PGI_BLASTOUT.txt", sep = "\t", header = FALSE)
summary(PGI_BLASTOUT$V5)
summary(PGI_BLASTOUT$V7)
hist(PGI_BLASTOUT$V5)
hist(PGI_BLASTOUT$V7, breaks = 100, xlim = c(0, 2000))

IBR31_BLASTOUT <- read.delim("IBR31_BLASTOUT.txt", sep = "\t", header = FALSE)
summary(IBR31_BLASTOUT$V5)
summary(IBR31_BLASTOUT$V7)
hist(IBR31_BLASTOUT$V5, breaks = 100)
hist(IBR31_BLASTOUT$V7, xlim = c(0,2500), breaks = 100)
lowQualityIBR31 <- subset(IBR31_BLASTOUT, IBR31_BLASTOUT$V5 < 200)
View(lowQualityIBR31)

IBR32_BLASTOUT <- read.delim("IBR32_BLASTOUT.txt", sep = "\t", header = FALSE)
summary(IBR32_BLASTOUT$V5)
summary(IBR32_BLASTOUT$V7)
hist(IBR32_BLASTOUT$V5, breaks = 100)
hist(IBR32_BLASTOUT$V7, breaks = 100)

LFY_barcodeBLASTOUT <- read.delim("LFY_barcodeBlastOut.txt", sep = "\t", header = FALSE)
summary(LFY_barcodeBLASTOUT$V5)
matches <- subset(LFY_barcodeBLASTOUT, LFY_barcodeBLASTOUT$V5 > 75)

PGI_primerBLASTOUT <- read.delim("PGI_barcodePrimerBlastOut.txt", sep = "\t", header = FALSE)
View(PGI_primerBLASTOUT)
hist(PGI_primerBLASTOUT$V7, breaks = 10)


### Flow Cytometry ###
flowcyt <- read.csv("PloidyCvalues.csv", header = TRUE)
setwd("~/Desktop/Isoetes_sandbox")

hist(flowcyt$C.value, breaks = 50, xlab = "C-value", main = "", prob = TRUE)
lines(density(flowcyt$C.value), lwd = 2)
diploids <- subset(flowcyt, flowcyt$Ploidy.Level == 22)
tetraploids <- subset(flowcyt, flowcyt$Ploidy.Level == 44)
hexaploids <- subset(flowcyt, flowcyt$Ploidy.Level == 66)
octoploids <- subset(flowcyt, flowcyt$Ploidy.Level == 88)
decaploids <- subset(flowcyt, flowcyt$Ploidy.Level == 110)
par(mfrow = c(1,1))
boxplot(diploids$C.value)
boxplot(tetraploids$C.value)
boxplot(diploids$C.value, tetraploids$C.value, hexaploids$C.value, octoploids$C.value, decaploids$C.value, ylab = "C-value", xlab = "Ploidy level including morpho-IDs", names = c("2x", "4x", "6x", "8x", "10x"), main =NA)

CvalVScopyNum <- lm(flowcyt$Ploidy.Level ~ flowcyt$C.value)
summary(CvalVScopyNum)

uniqueCval <- flowcyt[!duplicated(flowcyt[1:2]),]
UniqueCvalVScopyNum <- lm(uniqueCval$Ploidy.Level ~ uniqueCval$C.value)
summary(UniqueCvalVScopyNum)

Cval.pca <- prcomp(flowcyt, center = TRUE, scale. = TRUE)
summary(Cval.pca)

Cval.aov <- aov(flowcyt$C.value ~ flowcyt$Ploidy.Level)
summary(Cval.aov)
TukeyHSD(Cval.aov)
plot(Cval.aov, 1)



### PURC Testing Clusters ###
IBR31_clusters <- as.data.frame(read.csv("IBR3-1_clusters.csv", header = FALSE))
IBR32_clusters <- read.csv("IBR3-2_clusters.csv", header = FALSE)
PGI_clusters <- read.csv("PGI_clusters.csv", header = FALSE)
LFY_clusters <- read.csv("LFY_clusters.csv", header = FALSE)

PURC_counts <- as.data.frame(read.delim("PURC_counts.tsv", sep = "\t", header = TRUE))
PURC_counts[PURC_counts == 0] <- NA
boxplot(PURC_counts$IBR31, PURC_counts$IBR32, PURC_counts$PGIC, PURC_counts$LFY, names = c("IBR3-1", "IBR3-2", "pgiC", "LEAFY"), xlab = "Locus", ylab = "Num. of CCS Reads", main = "PURC Output")

hist(PURC_counts$IBR31, breaks = 100, prob = TRUE)
lines(density(PURC_counts$IBR31, na.rm = TRUE), lwd = 2)

hist(PURC_counts$IBR32, breaks = 100, prob = TRUE)
lines(density(PURC_counts$IBR32, na.rm = TRUE), lwd = 2)

hist(PURC_counts$LFY, breaks = 100, prob = TRUE)
lines(density(PURC_counts$LFY, na.rm = TRUE), lwd = 2)

hist(PURC_counts$PGIC, breaks = 100, prob = TRUE)
lines(density(PURC_counts$PGIC, na.rm = TRUE), lwd = 2)


boxplot(LFY_clusters$V1, PGI_clusters$V1, IBR31_clusters$V1, IBR32_clusters$V1)
IBR31_twoStdDev <- ggplot(data = IBR31_clusters, aes(x="", y=IBR31_clusters$V1)) + geom_boxplot()
meanIBR31 <- mean(IBR31_clusters$V1)
sigma1IBR31 <- meanIBR31 - sd(IBR31_clusters$V1)
sigma2IBR31 <- meanIBR31 + sd(IBR31_clusters$V1)
TwoSigma1IBR31 <- meanIBR31 - (2*sd(IBR31_clusters$V1))
TwoSigma2IBR31 <- meanIBR31 + (2*sd(IBR31_clusters$V1))
stdev <- sd(IBR31_clusters$V1)
IBR31_twoStdDev + stat_summary(fun.y = "meanIBR31", fun.ymin = "sigma1IBR31",  fun.ymax = "sigma2IBR31", geom = "bar")


### 2019 March 29 -- C-values CVs ###
setwd("~/Desktop/Isoetes_sandbox")
Cvals.subsamples <- read.csv("Cvalues_subsamples.csv", header = TRUE)
Cvals.popsamples <- read.csv("Cvalues_populationsamples.csv", header = TRUE)
Cvals.popsamples.rmOutliers <- read.csv("Cvalues_populationsamples_rmOutliers.csv", header = TRUE)
boxplot(Cvals.subsamples$CV, Cvals.popsamples$CV, names = c("Indv. Repl.", "Pop. Repl."), ylab = "CV", main = "C-Value Coefficients of Variation ")
boxplot(Cvals.subsamples$CV, Cvals.popsamples.rmOutliers$CV, names = c("Indv. Repl.", "Pop. Repl."))

t.test(x = Cvals.subsamples$CV, y = Cvals.popsamples.rmOutliers$CV)
plot(Cvals.popsamples$Replication~Cvals.popsamples$CV)
plot(Cvals.subsamples$Replication~Cvals.subsamples$CV)
plot(Cvals.popsamples.rmOutliers$CV~Cvals.popsamples.rmOutliers$Replication)
temp <- lm(Cvals.popsamples.rmOutliers$CV~as.numeric(Cvals.popsamples.rmOutliers$Replication))
summary(temp)

engelmannii <- read.csv("engelmannii.csv", header = TRUE)
mean(engelmannii$CV)
sd(engelmannii$MEAN)/mean(engelmannii$MEAN)*100
View(engelmannii)

