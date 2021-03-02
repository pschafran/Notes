setwd("/Users/Peter/Desktop/Isoetes_sandbox/assignClustersToDiploids")
polyploidSampleDistanceMatrix <- as.matrix(read.csv("Polyploid_Sample_Distance_PacBioPatristicDistances.csv", sep = ",", header = TRUE, row.names =1))
View(polyploidSampleDistanceMatrix)

library(phangorn)
njPolyploid <- NJ(polyploidSampleDistanceMatrix)
plot(njPolyploid, type = "phylo", cex = 0.25)
write.tree(njPolyploid, file = "NJ_polyploid_samples.tre")

upgmaPolyploid <- upgma(polyploidSampleDistanceMatrix)
plot(upgmaPolyploid, type = "u")
write.tree(upgmaPolyploid, file = "UPGMA_polyploid_samples.tre")

comparePhylo(njPolyploid, upgmaPolyploid, plot = TRUE, force.rooted = FALSE, use.edge.length = FALSE)

### 2019 Oct 28 ###
setwd("/Volumes/Samsung_T5/Isoetes_sandbox/")
library(phangorn)
polyploidSampleDistanceMatrix <- as.matrix(read.csv("Polyploid_Sample_Distance_pgiC_distancematrix.csv", sep = ",", header = TRUE, row.names =1))
dim(polyploidSampleDistanceMatrix)
pgiC_njPolyploid <- NJ(polyploidSampleDistanceMatrix)
plot(pgiC_njPolyploid, type = "phylo", cex = 0.25, label.offset = 0.001)
write.tree(njPolyploid, file = "NJ_pgiC_polyploid_samples.tre")

# IBR3-1
polyploidSampleDistanceMatrix <- as.matrix(read.csv("Polyploid_Sample_Distance_IBR3-1_distancematrix.csv", sep = ",", header = TRUE, row.names =1))
dim(polyploidSampleDistanceMatrix)
IBR31_njPolyploid <- NJ(polyploidSampleDistanceMatrix)
plot.phylo(IBR31_njPolyploid, type = "phylo", cex = 0.25, label.offset = 0.001, direction = "leftwards" )
write.tree(njPolyploid, file = "NJ_IBR3-1_polyploid_samples.tre")

# IBR3-2
polyploidSampleDistanceMatrix <- as.matrix(read.csv("Polyploid_Sample_Distance_IBR3-2_distancematrix.csv", sep = ",", header = TRUE, row.names =1))
dim(polyploidSampleDistanceMatrix)
IBR32_njPolyploid <- NJ(polyploidSampleDistanceMatrix)
upgmaPolyploid <- upgma(polyploidSampleDistanceMatrix)
plot(IBR32_njPolyploid, type = "phylo", cex = 0.25, label.offset = 0.001)
plot(upgmaPolyploid, type = "phylo", cex = 0.25)
write.tree(njPolyploid, file = "NJ_IBR3-2_polyploid_samples.tre")
comparePhylo(upgmaPolyploid, njPolyploid, plot = TRUE)

# LEAFY

polyploidSampleDistanceMatrix <- as.matrix(read.csv("/Volumes/Samsung_T5/IsoetesDNA/PacBio/assignClustersToDiploids/LEAFY/Polyploid_Sample_Distance_PacBioPatristicDistances.csv", sep = ",", header = TRUE, row.names =1))
dim(polyploidSampleDistanceMatrix)
LFY_njPolyploid <- NJ(polyploidSampleDistanceMatrix)
upgmaPolyploid <- upgma(polyploidSampleDistanceMatrix)
plot(LFY_njPolyploid, type = "phylo", cex = 0.25, label.offset = 0.001)
plot(upgmaPolyploid, type = "phylo", cex = 0.25)
write.tree(LFY_njPolyploid, file = "NJ_LEAFY_polyploid_samples.tre")
comparePhylo(upgmaPolyploid, njPolyploid, plot = TRUE)


# Create lists of just sample IDs
for (tip in LFY_njPolyploid$tip.label){
  if (exists("LFY_samplelist")){
    LFY_samplelist <- c(LFY_samplelist, strsplit(tip, "_")[[1]][3])
  } else (LFY_samplelist <- strsplit(tip, "_")[[1]][3])
}

for (tip in pgiC_njPolyploid$tip.label){
  if (exists("pgiC_samplelist")){
    pgiC_samplelist <- c(pgiC_samplelist, strsplit(tip, "_")[[1]][3])
  } else (pgiC_samplelist <- strsplit(tip, "_")[[1]][3])
}

for (tip in IBR31_njPolyploid$tip.label){
  if (exists("IBR31_samplelist")){
    IBR31_samplelist <- c(IBR31_samplelist, strsplit(tip, "_")[[1]][3])
  } else (IBR31_samplelist <- strsplit(tip, "_")[[1]][3])
}

for (tip in IBR32_njPolyploid$tip.label){
  if (exists("IBR32_samplelist")){
    IBR32_samplelist <- c(IBR32_samplelist, strsplit(tip, "_")[[1]][3])
  } else (IBR32_samplelist <- strsplit(tip, "_")[[1]][3])
}

LFY_samplelist
pgiC_samplelist
IBR31_samplelist
IBR32_samplelist


# Syncronize samples across trees

for (tip in pgiC_njPolyploid$tip.label){
  splitTip <- strsplit(tip, "_")[[1]][3]
  if ((splitTip %in% IBR31_samplelist) & (splitTip %in% IBR32_samplelist) & (splitTip %in% LFY_samplelist)) {
  } else {
    pgiC_njPolyploid <- drop.tip(pgiC_njPolyploid, tip)
  }
  print(length(pgiC_njPolyploid$tip.label))
}

for (tip in IBR31_njPolyploid$tip.label){
  splitTip <- strsplit(tip, "_")[[1]][3]
  if ((splitTip %in% LFY_samplelist) & (splitTip %in% IBR32_samplelist) & (splitTip %in% pgiC_samplelist)) {
  } else {
    IBR31_njPolyploid <- drop.tip(IBR31_njPolyploid, tip)
  }
  print(length(IBR31_njPolyploid$tip.label))
}

for (tip in IBR32_njPolyploid$tip.label){
  splitTip <- strsplit(tip, "_")[[1]][3]
  if ((splitTip %in% IBR31_samplelist) & (splitTip %in% pgiC_samplelist) & (splitTip %in% LFY_samplelist)) {
  } else {
    IBR32_njPolyploid <- drop.tip(IBR32_njPolyploid, tip)
  }
  print(length(IBR32_njPolyploid$tip.label))
}

for (tip in LFY_njPolyploid$tip.label){
  splitTip <- strsplit(tip, "_")[[1]][3]
  if ((splitTip %in% IBR31_samplelist) & (splitTip %in% pgiC_samplelist) & (splitTip %in% IBR32_samplelist)) {
  } else {
    LFY_njPolyploid <- drop.tip(LFY_njPolyploid, tip)
  }
  print(length(LFY_njPolyploid$tip.label))
}

# Syncronize tips only IBR and pgiC

for (tip in pgiC_njPolyploid$tip.label){
  if ((tip %in% IBR31_njPolyploid$tip.label) & (tip %in% IBR32_njPolyploid$tip.label)) {
  } else {
    pgiC_njPolyploid <- drop.tip(pgiC_njPolyploid, tip)
  }
  print(length(pgiC_njPolyploid$tip.label))
}

for (tip in IBR31_njPolyploid$tip.label){
  if ((tip %in% IBR32_njPolyploid$tip.label) & (tip %in% pgiC_njPolyploid$tip.label)) {
  } else {
    IBR31_njPolyploid <- drop.tip(IBR31_njPolyploid, tip)
  }
  print(length(IBR31_njPolyploid$tip.label))
}

for (tip in IBR32_njPolyploid$tip.label){
  if ((tip %in% IBR31_njPolyploid$tip.label) & (tip %in% pgiC_njPolyploid$tip.label)) {
  } else {
    IBR32_njPolyploid <- drop.tip(IBR32_njPolyploid, tip)
  }
  print(length(IBR32_njPolyploid$tip.label))
}



length(pgiC_njPolyploid$tip.label)
length(IBR31_njPolyploid$tip.label)
length(IBR32_njPolyploid$tip.label)
length(LFY_njPolyploid$tip.label)

write.tree(pgiC_njPolyploid, file = "NJ_pgiC_reduced_155_samples.tre")
write.tree(IBR31_njPolyploid, file = "NJ_IBR31_reduced_155_samples.tre")
write.tree(IBR32_njPolyploid, file = "NJ_IBR32_reduced_155_samples.tre")
write.tree(LFY_njPolyploid, file = "NJ_LEAFY_reduced_polyploid_samples.tre")
LFY_njPolyploid <- read.tree("NJ_LEAFY_reduced_polyploid_samples.tre")

pgiC_njPolyploid$tip.label
IBR31_njPolyploid$tip.label
IBR32_njPolyploid$tip.label
LFY_njPolyploid$tip.label

setdiff(pgiC_njPolyploid$tip.label, IBR31_njPolyploid$tip.label)
setdiff(IBR32_njPolyploid$tip.label, IBR31_njPolyploid$tip.label)
setdiff(LFY_njPolyploid$tip.label, pgiC_njPolyploid$tip.label)
setdiff(pgiC_njPolyploid$tip.label, LFY_njPolyploid$tip.label)
comparePhylo(LFY_njPolyploid, IBR32_njPolyploid)
IBR32_njPolyploid <- drop.tip(IBR32_njPolyploid, "Isoetes_hyemalis_Schafran120_1")

consensusTree <- consensus(pgiC_njPolyploid, IBR31_njPolyploid, IBR32_njPolyploid, LFY_njPolyploid, p = 0.5, check.labels = TRUE)
plot.phylo(consensusTree, "u", lab4ut = "axial", cex = 0.5)

trees <- c(pgiC_njPolyploid, IBR31_njPolyploid, IBR32_njPolyploid, LFY_njPolyploid)
consensusNet <- consensusNet(trees, p = 0.3)
plot(consensusNet, "2D", cex = 0.5)

consensusTree3loci <- consensus(pgiC_njPolyploid, IBR31_njPolyploid, IBR32_njPolyploid, p = 0.5, check.labels = TRUE)
plot(consensusTree3loci, "u", lab4ut = "axial", cex = 0.5)
write.tree(consensusTree3loci, file = "PacBio_Consensus_3loci.tre")

threetrees <- c(pgiC_njPolyploid, IBR31_njPolyploid, IBR32_njPolyploid) 
consensusNet3loci <- consensusNet(threetrees, p = 0.5)
plot(consensusNet3loci, "2D", cex = 0.5)
