
require(ape)
require(phangorn)
require(phytools)

setwd("/Volumes/Samsung_T5/IsoetesDNA/GoFlag/HybPiper/")

hybPiper_RawConcat <- read.tree("./noSinensisRef/HybPiper_SingleCopy_ConcatTree.newick")
hybPiper_RawAstral <- read.tree("./noSinensisRef/HybPiper_SingleCopy_ASTRAL.tre")

plot(hybPiper_RawAstral)
plot(hybPiper_RawConcat)
plot.phylo(hybPiper_RawAstral, use.edge.length = FALSE)
plot.phylo(hybPiper_RawConcat, use.edge.length = FALSE)


hybPiper_RawAstral$tip.label

par(mfrow = c(1,1))
hybPiper_noFilter <- comparePhylo(hybPiper_RawAstral, hybPiper_RawConcat, plot = TRUE, force.rooted = FALSE, use.edge.length = FALSE)
dist.topo(hybPiper_RawAstral, hybPiper_RawConcat, method = "PH85")

edge.lab <- createLabel(hybPiper_RawAstral, hybPiper_RawConcat, hybPiper_RawConcat$edge[,2], "edge")
edge.lab
plot(hybPiper_RawConcat, "u", main = "Concatenated Unfiltered", cex = 0.25)
edgelabels(hybPiper_RawConcat$edge[,2], col = "blue", cex = 0.5, frame = "none")
plot(hybPiper_RawAstral, "u", main = "ASTRAL Unfiltered", cex = 0.5)
edgelabels(hybPiper_RawAstral$edge[,2], col = "blue", cex = 0.5, frame = "none")
edgelabels(hybPiper_RawAstral, col = "red", cex = 1, frame = "rect")
edge.col <- rep("black", nrow(hybPiper_RawAstral$edge))
edge.col[ is.na(edge.lab) ] <- "red"
x <- plot(hybPiper_RawAstral, "phylogram", edge.color = edge.col, cex=.5, main = "Unshared Branches", use.edge.length = FALSE)

### HybPiper Diploids Target Region
setwd("/Volumes/Samsung_T5/IsoetesDNA/GoFlag/HybPiper/diploidsOnly/targetRegions/")

## IQ-TREE
hybPiper_target_strict_100pct_contrees <- read.tree("HybPiper.Diploid.TargetRegion.STRICT.100pct.IQTREE.contrees.ASTRAL")
hybPiper_target_strict_100pct_trees <- read.tree("HybPiper.Diploid.TargetRegion.STRICT.100pct.IQTREE.trees.ASTRAL")
hybPiper_target_strict_100pct_contrees <- unroot(hybPiper_target_strict_100pct_contrees)
hybPiper_target_strict_100pct_trees <- unroot(hybPiper_target_strict_100pct_trees)

# Compare consensus trees vs. likelihod trees
plot(hybPiper_target_strict_100pct_trees)
plot(hybPiper_target_strict_100pct_contrees)
comparePhylo(hybPiper_target_strict_100pct_contrees, hybPiper_target_strict_100pct_trees, plot = TRUE, force.rooted = FALSE, use.edge.length = FALSE)
dist.topo(hybPiper_target_strict_100pct_contrees, hybPiper_target_strict_100pct_trees, method = "PH85")

## RAxML
HybPiper.Diploid.TargetRegion.AUTO.RAxML.ASTRAL	<- read.tree("HybPiper.Diploid.TargetRegion.AUTO.RAxML.ASTRAL")
HybPiper.Diploid.TargetRegion.NONE.RAxML.ASTRAL <- read.tree("HybPiper.Diploid.TargetRegion.NONE.RAxML.ASTRAL")
HybPiper.Diploid.TargetRegion.STRICT.RAxML.ASTRAL <- read.tree("HybPiper.Diploid.TargetRegion.STRICT.RAxML.ASTRAL")
HybPiper.Diploid.TargetRegion.SUPERSTRICT.RAxML.ASTRAL <- read.tree("HybPiper.Diploid.TargetRegion.SUPERSTRICT.RAxML.ASTRAL")

dist.topo(unroot(HybPiper.Diploid.TargetRegion.NONE.RAxML.ASTRAL), unroot(HybPiper.Diploid.TargetRegion.AUTO.RAxML.ASTRAL), method = "PH85")
dist.topo(unroot(HybPiper.Diploid.TargetRegion.NONE.RAxML.ASTRAL), unroot(HybPiper.Diploid.TargetRegion.STRICT.RAxML.ASTRAL), method = "PH85")
dist.topo(unroot(HybPiper.Diploid.TargetRegion.NONE.RAxML.ASTRAL), unroot(HybPiper.Diploid.TargetRegion.SUPERSTRICT.RAxML.ASTRAL), method = "PH85")
dist.topo(unroot(HybPiper.Diploid.TargetRegion.AUTO.RAxML.ASTRAL), unroot(HybPiper.Diploid.TargetRegion.STRICT.RAxML.ASTRAL), method = "PH85")
dist.topo(unroot(HybPiper.Diploid.TargetRegion.AUTO.RAxML.ASTRAL), unroot(HybPiper.Diploid.TargetRegion.SUPERSTRICT.RAxML.ASTRAL), method = "PH85")
dist.topo(unroot(HybPiper.Diploid.TargetRegion.STRICT.RAxML.ASTRAL), unroot(HybPiper.Diploid.TargetRegion.SUPERSTRICT.RAxML.ASTRAL), method = "PH85")
comparePhylo(HybPiper.Diploid.TargetRegion.STRICT.RAxML.ASTRAL, HybPiper.Diploid.TargetRegion.SUPERSTRICT.RAxML.ASTRAL, plot = TRUE, force.rooted = FALSE, use.edge.length = FALSE)

ultrmtrc.HybPiper.Diploid.Target.SUPER.RAxML.ASTRAL <- (force.ultrametric(HybPiper.Diploid.TargetRegion.SUPERSTRICT.RAxML.ASTRAL, method = "extend"))
ultrmtrc.HybPiper.Diploid.Target.STRICT.RAxML.ASTRAL <- (force.ultrametric(HybPiper.Diploid.TargetRegion.STRICT.RAxML.ASTRAL, method = "extend"))

plot(ultrmtrc.HybPiper.Diploid.Target.STRICT.RAxML.ASTRAL, "u", cex = 0.5, main = "Strict Alignment Trimming")
plot(ultrmtrc.HybPiper.Diploid.Target.SUPER.RAxML.ASTRAL, "u", main = "Superstrict Alignment Trimming", cex = 0.5)
edge.lab <- createLabel(ultrmtrc.HybPiper.Diploid.Target.STRICT.RAxML.ASTRAL, ultrmtrc.HybPiper.Diploid.Target.SUPER.RAxML.ASTRAL, ultrmtrc.HybPiper.Diploid.Target.SUPER.RAxML.ASTRAL$edge[,2], "edge")
edge.col <- rep("black", nrow(ultrmtrc.HybPiper.Diploid.Target.STRICT.RAxML.ASTRAL$edge))
edge.col[ is.na(edge.lab) ] <- "white"
x <- plot(ultrmtrc.HybPiper.Diploid.Target.SUPER.RAxML.ASTRAL, "u", edge.color = edge.col, cex=0.75, main = "Strict vs. Superstrict Unshared Branches", use.edge.length = FALSE)

edge.lab <- createLabel(HybPiper.Diploid.TargetRegion.NONE.RAxML.ASTRAL, HybPiper.Diploid.TargetRegion.SUPERSTRICT.RAxML.ASTRAL, HybPiper.Diploid.TargetRegion.SUPERSTRICT.RAxML.ASTRAL$edge[,2], "edge")
edge.col <- rep("black", nrow(HybPiper.Diploid.TargetRegion.NONE.RAxML.ASTRAL$edge))
edge.col[ is.na(edge.lab) ] <- "white"
x <- plot(HybPiper.Diploid.TargetRegion.SUPERSTRICT.RAxML.ASTRAL, "u", edge.color = edge.col, cex=0.75, main = "No Filter vs. Superstrict Unshared Branches", use.edge.length = FALSE)

edge.lab <- createLabel(HybPiper.Diploid.TargetRegion.AUTO.RAxML.ASTRAL, HybPiper.Diploid.TargetRegion.STRICT.RAxML.ASTRAL, HybPiper.Diploid.TargetRegion.STRICT.RAxML.ASTRAL$edge[,2], "edge")
edge.col <- rep("black", nrow(HybPiper.Diploid.TargetRegion.AUTO.RAxML.ASTRAL$edge))
edge.col[ is.na(edge.lab) ] <- "white"
x <- plot(HybPiper.Diploid.TargetRegion.STRICT.RAxML.ASTRAL, "u", edge.color = edge.col, cex=0.75, main = "Auto vs. Strict Unshared Branches", use.edge.length = FALSE)

edge.lab <- createLabel(HybPiper.Diploid.TargetRegion.NONE.RAxML.ASTRAL, HybPiper.Diploid.TargetRegion.STRICT.RAxML.ASTRAL, HybPiper.Diploid.TargetRegion.STRICT.RAxML.ASTRAL$edge[,2], "edge")
edge.col <- rep("black", nrow(HybPiper.Diploid.TargetRegion.NONE.RAxML.ASTRAL$edge))
edge.col[ is.na(edge.lab) ] <- "white"
x <- plot(HybPiper.Diploid.TargetRegion.STRICT.RAxML.ASTRAL, "u", edge.color = edge.col, cex=0.75, main = "No Filter vs. Strict Unshared Branches", use.edge.length = FALSE)

# RAxML Concat Trees

HybPiper.Diploid.TargetRegion.NONE.RAxML.concat <- read.tree("HybPiper.Diploid.TargetRegion.NONE.RAxML.concat")
HybPiper.Diploid.TargetRegion.AUTO.RAxML.concat <- read.tree("HybPiper.Diploid.TargetRegion.AUTO.RAxML.concat")
HybPiper.Diploid.TargetRegion.STRICT.RAxML.concat <- read.tree("HybPiper.Diploid.TargetRegion.STRICT.RAxML.concat")
HybPiper.Diploid.TargetRegion.SUPERSTRICT.RAxML.concat <- read.tree("HybPiper.Diploid.TargetRegion.SUPERSTRICT.RAxML.concat")

dist.topo(unroot(HybPiper.Diploid.TargetRegion.NONE.RAxML.concat), unroot(HybPiper.Diploid.TargetRegion.AUTO.RAxML.concat), method = "PH85")
dist.topo(unroot(HybPiper.Diploid.TargetRegion.NONE.RAxML.concat), unroot(HybPiper.Diploid.TargetRegion.STRICT.RAxML.concat), method = "PH85")
dist.topo(unroot(HybPiper.Diploid.TargetRegion.NONE.RAxML.concat), unroot(HybPiper.Diploid.TargetRegion.SUPERSTRICT.RAxML.concat), method = "PH85")
dist.topo(unroot(HybPiper.Diploid.TargetRegion.AUTO.RAxML.concat), unroot(HybPiper.Diploid.TargetRegion.STRICT.RAxML.concat), method = "PH85")
dist.topo(unroot(HybPiper.Diploid.TargetRegion.AUTO.RAxML.concat), unroot(HybPiper.Diploid.TargetRegion.SUPERSTRICT.RAxML.concat), method = "PH85")
dist.topo(unroot(HybPiper.Diploid.TargetRegion.STRICT.RAxML.concat), unroot(HybPiper.Diploid.TargetRegion.SUPERSTRICT.RAxML.concat), method = "PH85")

par(mfrow=c(1,1))



### HybPiper Diploids Supercontigs
setwd("/Volumes/Samsung_T5/IsoetesDNA/GoFlag/HybPiper/diploidsOnly/supercontigs/50pct")
hybPiper_supercontigs_auto_50pct_trees <- read.tree("HybPiper.Diploids.Supercontigs.AUTO.50pct.trees.ASTRAL")
plot(hybPiper_supercontigs_auto_50pct_trees)


### HybPiper Diploids Introns
setwd("/Volumes/Samsung_T5/IsoetesDNA/GoFlag/HybPiper/diploidsOnly/introns/100pct")
hybPiper_introns_100pct_trees <- read.tree("HybPiper.Diploid.Introns.100pct.IQTREE.trees.ASTRAL")
hybPiper_introns_100pct_contrees <- read.tree("HybPiper.Diploid.Introns.100pct.IQTREE.contrees.ASTRAL")
plot(hybPiper_introns_100pct_trees, "c")
plot(hybPiper_introns_100pct_contrees, "p")
comparePhylo(hybPiper_introns_100pct_contrees, hybPiper_introns_100pct_trees, plot = TRUE, force.rooted = FALSE, use.edge.length = FALSE)


### Pairwise Comparisons b/w ASTRAL IQTREE target region results 

# 50% missing data rows




### Target Recovery Phyluce vs HybPiper vs HybPhyloMaker

setwd("/Volumes/Samsung_T5/IsoetesDNA/GoFlag/IsoetesPipeline")
phyluce.percRecovered <- read.table("Phyluce-diploids-percRecovered.txt", row.names = 1)
hybpiper.percRecovered <- read.table("HybPiper-diploids-percRecovered.txt", row.names = 1)
hybphylo.percRecovered <- read.table("HybPhyloMaker-diploids-percRecovered.txt", row.names = 1)
goflag.percRecovered <- read.table("GoFlag_percRecovered.txt", row.names = 1)

length(phyluce.percRecovered[,1])
length(hybpiper.percRecovered[,1])
length(hybphylo.percRecovered[,1])

merged.percRecovered <- merge(phyluce.percRecovered, hybpiper.percRecovered, by="row.names", all=TRUE)
merged.percRecovered[is.na(merged.percRecovered)] <- 0
rownames(merged.percRecovered) <- merged.percRecovered[,1]
merged.percRecovered <- subset(merged.percRecovered, select = c(V2.x,V2.y))

merged.percRecovered <- merge(merged.percRecovered, hybphylo.percRecovered, by="row.names", all=TRUE)
merged.percRecovered[is.na(merged.percRecovered)] <- 0
rownames(merged.percRecovered) <- merged.percRecovered[,1]
merged.percRecovered <- subset(merged.percRecovered, select = c(V2.x,V2.y,V2))

merged.percRecovered <- merge(merged.percRecovered, goflag.percRecovered, by="row.names", all=TRUE)
merged.percRecovered[is.na(merged.percRecovered)] <- 0
rownames(merged.percRecovered) <- merged.percRecovered[,1]
colnames(merged.percRecovered) <- c("rowNames", "Phyluce", "HybPiper", "HybPhyloMaker","GoFlag")
merged.percRecovered <- subset(merged.percRecovered, select = c("Phyluce", "HybPiper", "HybPhyloMaker","GoFlag"))


View(merged.percRecovered)
heatmap.2(as.matrix(merged.percRecovered), Rowv = FALSE, Colv = FALSE, dendrogram= c('none'), 
          density.info = 'none', trace = 'none', 
          cexCol = 0.8,
          lmat=rbind( 3:4,2:1 ), lhei=c(1,4), lwid=c(1,4), 
          col = gray.colors(10,start=0.9,end=0.3,gamma=2.2),
          key.title = NA, key.xlab = "Percent Recovery", keysize = 1
          )
require(lattice)
rgbPalette <- colorRampPalette(c("red","green","blue","purple"), space = "Lab")
levelplot(as.matrix(merged.percRecovered), scale=list(x=list(rot=45)), aspect = "fill", xlab = "", ylab = "", col.regions = rgbPalette, main = "Percent Recovery of Targets", cex = 0.5)
colnames(goflag.percRecovered$V2) < "GoFlag"
levelplot(as.matrix(goflag.percRecovered), scale=list(x=list(rot=45)), aspect = "fill", xlab = "", ylab = "", col.regions = rgbPalette, main = "GoFlag Percent Recovery of Targets", cex = 0.5)
boxplot(goflag.percRecovered$V2, main = "GoFlag Percent Recovery of Targets", ylab = "Percentage of Targets Recovered")



### HybPiper Trees
require(phytools)
require(phangorn)
require(ks)
setwd("/Volumes/Samsung_T5/IsoetesDNA/GoFlag/IsoetesPipeline/HybPiper")
hybpipertrees <- read.tree(file.choose())
for (i in 1:9) {
hybpipertrees[[i]]$edge.length[is.na(hybpipertrees[[i]]$edge.length)] <- 0
}
qpd.tree <- averageTree(hybpipertrees, method = "quadratic.path.difference")
tmp<-list(qpd.tree)
class(tmp)<-"multiPhylo"
obj<-c(hybpipertrees,tmp)
D<-phytools:::qpd(obj,obj)
MDS<-cmdscale(D)
plot(MDS,xlab="MDS axis 1",ylab="MDS axis 2")
for(i in 1:nrow(MDS)) 
  lines(c(MDS[i,1],MDS[nrow(MDS),1]),c(MDS[i,2],MDS[nrow(MDS),2]),
        col="grey")
points(MDS[1,1],MDS[1,2],pch=21,bg="blue", cex = 1)
points(MDS[2,1],MDS[2,2],pch=22,bg="blue", cex = 1)
points(MDS[3,1],MDS[3,2],pch=23,bg="blue", cex = 1)
points(MDS[4,1],MDS[4,2],pch=21,bg="green", cex = 1)
points(MDS[5,1],MDS[5,2],pch=22,bg="green", cex = 1)
points(MDS[6,1],MDS[6,2],pch=23,bg="green", cex = 1)
points(MDS[7,1],MDS[7,2],pch=21,bg="red", cex = 1)
points(MDS[8,1],MDS[8,2],pch=22,bg="red", cex = 1)
points(MDS[9,1],MDS[9,2],pch=23,bg="red", cex = 1)
points(MDS[10,1],MDS[10,2],pch=21,bg="yellow", cex = 1)
points(MDS[11,1],MDS[11,2],pch=22,bg="yellow", cex = 1)
points(MDS[12,1],MDS[12,2],pch=23,bg="yellow", cex = 1)
points(x=MDS[nrow(MDS),1],y=MDS[nrow(MDS),2],pch=21,bg="black")
abline(v=0,h=0,lty="dotted")
legend(0.5,0.1, legend = c("MAFFT", "Auto", "Strict", "Superstrict"), col = c("green", "blue", "red", "gold"), pch = 19, box.lty = 0)
legend(1.5,0.1, legend = c("100% data", "75% data", "50% data"), pch = c(21, 22, 23), box.lty = 0)

d<-cmdscale(D,k=1)[,1]
h<-sapply(hybpipertrees,function(x) max(nodeHeights(x)))
pd<-density(d,bw=1)
plot(pd$x,pd$y,lwd=2,type="l",ylab="density",
     xlab="QPD distance (one-dimensional MDS)")

pd<-density(d)
plot(pd$x,pd$y,lwd=2,type="l",ylab="density",
     xlab="QPD distance (one-dimensional MDS)")
lines(rep(d[121],2),par()$usr[3:4],lty="dashed",col="red")
arrows(x0=d[121],y0=par()$usr[4],x1=d[121],
       y1=pd$y[which(abs(pd$x-d[121])==min(abs(pd$x-d[121])))],
       lwd=2,col="red",length=0.15,angle=20)
text(x=d[121],y=0.98*par()$usr[4],"average tree",pos=4,cex=0.8)

MDS<-cmdscale(D)
dim(MDS)
pd<-kde(MDS)
plot(pd,xlab="canonical axis 1",ylab="canonical axis 2")
points(MDS[1:39],pch=19)
points(MDS[39,1],MDS[39,2],pch=21,cex=1.25,bg="grey")
rect(MDS[39,1]+0.5*strwidth("W"),
     MDS[39,2]-0.5*0.9*strheight("W"),
     MDS[39,1]+0.5*strwidth("W")+0.9*strwidth("average tree"),
     MDS[39,2]+0.5*0.9*strheight("W"),border="transparent",
     col="white")
text(MDS[39,1],MDS[39,2],"average tree",pos=4,cex=0.8)


nucl.Tree <- read.tree(file.choose())
plast.Tree <- read.tree(file.choose())
unroot(nucl.Tree)
unroot(plast.Tree)
comparePhylo(nucl.Tree,plast.Tree, plot = TRUE)
plot(obj<-cophylo(nucl.Tree,plast.Tree), cex = 0.5)
plot(consensus(nucl.Tree, plast.Tree), "u")

all.Astral <- read.tree(file.choose())
unroot(all.Astral)
nucl.conTree <- consensus(all.Astral, p = 0.5)
plot(nucl.conTree, "u")
unroot(nucl.conTree)
cophyloplot(nucl.conTree,plast.Tree)
plot(obj<-cophylo(nucl.conTree,plast.Tree), cex = 0.5)
comparePhylo(nucl.conTree,plast.Tree, plot = TRUE, use.edge.length = FALSE)
par(mfrow = c(1,2))
plot(nucl.conTree, "u", cex = 0.5, use.edge.length = FALSE)
plot(plast.Tree, "u", cex = 0.5, use.edge.length = FALSE, rotate.tree = 45)


superstrict.Astral <- read.tree(file.choose())
superstrict.Concat <- read.tree(file.choose())
unroot(superstrict.Astral)
unroot(superstrict.Concat)
par(mfrow = c(1,1))
plot(superstrict.Astral, "u", cex = 0.5, use.edge.length = FALSE, rotate.tree = -70)
plot(plast.Tree, "p", cex = 0.5, use.edge.length = TRUE, rotate.tree = 45)
plot(superstrict.Concat, "u", cex = 0.5, use.edge.length = TRUE, rotate.tree = 90)



### Phased Tetraploids

strict.phased.trees <- read.tree(file.choose())

completeTrees <- list()
j <- 0
for (i in 1:174) {
  strict.phased.trees[[i]]$edge.length[is.na(strict.phased.trees[[i]]$edge.length)] <- 0
  if (length(strict.phased.trees[[i]]$tip.label) > 18) {
    j <- j+1
    print(length(strict.phased.trees[[i]]$tip.label))
    print(j)
    completeTrees[j] <- strict.phased.trees[i]
  }
}
qpd.tree <- averageTree(completeTrees, method = "quadratic.path.difference")
tmp<-list(qpd.tree)
class(tmp)<-"multiPhylo"
obj<-c(completeTrees,tmp)
D<-phytools:::qpd(obj,obj)
MDS<-cmdscale(D)
plot(MDS,xlab="MDS axis 1",ylab="MDS axis 2")
for(i in 1:nrow(MDS)) 
  lines(c(MDS[i,1],MDS[nrow(MDS),1]),c(MDS[i,2],MDS[nrow(MDS),2]),
        col="grey")
points(MDS[1:19,1],MDS[1:19,2],pch=21,bg="gray", cex = 1)
points(x=MDS[nrow(MDS),1],y=MDS[nrow(MDS),2],pch=21,bg="black")
abline(v=0,h=0,lty="dotted")

d<-cmdscale(D,k=1)[,1]
h<-sapply(completeTrees,function(x) max(nodeHeights(x)))
pd<-density(d,bw=1)
plot(pd$x,pd$y,lwd=2,type="l",ylab="density",
     xlab="QPD distance (one-dimensional MDS)")

pd<-density(d)
plot(pd$x,pd$y,lwd=2,type="l",ylab="density",
     xlab="QPD distance (one-dimensional MDS)")
lines(rep(d[121],2),par()$usr[3:4],lty="dashed",col="red")
arrows(x0=d[121],y0=par()$usr[4],x1=d[121],
       y1=pd$y[which(abs(pd$x-d[121])==min(abs(pd$x-d[121])))],
       lwd=2,col="red",length=0.15,angle=20)
text(x=d[121],y=0.98*par()$usr[4],"average tree",pos=4,cex=0.8)

MDS<-cmdscale(D)
dim(MDS)
pd<-kde(MDS)
plot(pd,xlab="canonical axis 1",ylab="canonical axis 2")
points(MDS[1:19],pch=19)
points(MDS[20,1],MDS[20,2],pch=21,cex=1.25,bg="grey")
rect(MDS[20,1]+0.5*strwidth("W"),
     MDS[20,2]-0.5*0.9*strheight("W"),
     MDS[20,1]+0.5*strwidth("W")+0.9*strwidth("average tree"),
     MDS[20,2]+0.5*0.9*strheight("W"),border="transparent",
     col="white")
text(MDS[20,1],MDS[20,2],"average tree",pos=4,cex=0.8)

plot(completeTrees[[1]])

fastDist<-function(tree,sp1,sp2){
  fastHeight(tree,sp1,sp1)+fastHeight(tree,sp2,sp2)-
    2*fastHeight(tree,sp1,sp2)
}

######################################################

superstrict.phased.trees <- read.tree(file.choose())
View(superstrict.phased.trees)
appTrees <- list()
j <- 0
for (i in 1:36) {
  x <- c("Isoetes-engelmannii_WA10_S1_L001", "I_appPhase0", "Isoetes-valida_WC11_L001", "I_appPhase1")
  superstrict.phased.trees[[i]]$edge.length[is.na(superstrict.phased.trees[[i]]$edge.length)] <- 0
  if (all(is.element(x,superstrict.phased.trees[[i]]$tip.label))) {
    j <- j+1
    print(i)
    print(is.element(x,superstrict.phased.trees[[i]]$tip.label))
    appTrees[j] <- superstrict.phased.trees[i]
  }
}
View(appTrees)
appMatrix <- matrix(data = NA, ncol = 8, nrow = 20)
colnames(appMatrix) <- c("eng-app0", "val-app0", "app0-app1", "val-app1", "eng-app1", "app0-app1", "App0 Min Distance", "App1 Min Distance")
j <- 1
for (i in seq(1,10,1)) {
appMatrix[j,1] <- fastDist(appTrees[[i]], "Isoetes-engelmannii_WA10_S1_L001", "I_appPhase0")
appMatrix[j,2] <- fastDist(appTrees[[i]], "Isoetes-valida_WC11_L001", "I_appPhase0")
appMatrix[j,3] <- fastDist(appTrees[[i]], "I_appPhase0", "I_appPhase1")
appMatrix[j+1,4] <- fastDist(appTrees[[i]], "Isoetes-valida_WC11_L001", "I_appPhase1")
appMatrix[j+1,5] <- fastDist(appTrees[[i]], "Isoetes-engelmannii_WA10_S1_L001", "I_appPhase1")
appMatrix[j+1,6] <- fastDist(appTrees[[i]], "I_appPhase0", "I_appPhase1")
appMatrix[j,7] <- colnames(appMatrix)[which.min(appMatrix[i,])]
appMatrix[j+1,8] <- colnames(appMatrix)[which.min(appMatrix[i+1,])]
j <- j + 2
}

View(appMatrix)

phaseSums <- matrix(data = NA, nrow = 1, ncol = 6)
colnames(phaseSums) <- c("eng-app0", "val-app0", "app0-app1", "eng-app1", "val-app1", "app0-app1")
phaseSums[1,1] <- sum(appMatrix[,7] == "eng-app0", na.rm = TRUE)
phaseSums[1,2] <- sum(appMatrix[,7] == "val-app0", na.rm = TRUE)
phaseSums[1,3] <- sum(appMatrix[,7] == "app0-app1", na.rm = TRUE)
phaseSums[1,4] <- sum(appMatrix[,8] == "eng-app1", na.rm = TRUE)
phaseSums[1,5] <- sum(appMatrix[,8] == "val-app1", na.rm = TRUE)
phaseSums[1,6] <- sum(appMatrix[,8] == "app0-app1", na.rm = TRUE)

View(phaseSums)
x <- barplot(phaseSums, xaxt = "n", ylab = "Number of genes supporting relationship")
labs <- colnames(phaseSums)
text(cex=1, x=x-.25, y=-0.5, labs, xpd=TRUE, srt=45)


fastDist(plast.Tree, "Isoetes-flaccida_WD11_L001", "Isoetes-melanopoda_WH11_L001")
fastDist(plast.Tree, "Isoetes-mattaponica_WB11_L001", "Isoetes-silvatica_WG11_L001")
fastDist(plast.Tree, "Isoetes-flaccida_WD11_L001", "Isoetes-valida_WC11_L001")
(fastDist(plast.Tree, "Isoetes-flaccida_WD11_L001", "Isoetes-melanopoda_WH11_L001")/mean(fastDist(plast.Tree, "Isoetes-mattaponica_WB11_L001", "Isoetes-silvatica_WG11_L001") + fastDist(plast.Tree, "Isoetes-flaccida_WD11_L001", "Isoetes-valida_WC11_L001") + fastDist(plast.Tree, "Isoetes-flaccida_WD11_L001", "Isoetes-melanopoda_WH11_L001")))*100

fastDist(superstrict.Concat, "Isoetes-flaccida_WD11_L001", "Isoetes-flaccida_WE11_L001")
fastDist(superstrict.Concat, "Isoetes-melanopoda_WH11_L001", "Isoetes-lithophila_WB12_L001")
fastDist(superstrict.Concat, "Isoetes-mattaponica_WB11_L001", "Isoetes-silvatica_WG11_L001")
(fastDist(superstrict.Concat, "Isoetes-flaccida_WD11_L001", "Isoetes-flaccida_WE11_L001")/mean(fastDist(superstrict.Concat, "Isoetes-melanopoda_WH11_L001", "Isoetes-lithophila_WB12_L001") + fastDist(superstrict.Concat, "Isoetes-melanopoda_WH11_L001", "Isoetes-lithophila_WB12_L001") + fastDist(superstrict.Concat, "Isoetes-mattaponica_WB11_L001", "Isoetes-silvatica_WG11_L001")))*100



### Adegenet Testing
setwd("/Volumes/Samsung_T5/IsoetesDNA/GoFlag/IsoetesPipeline/HybPiper/diploid-results/dna/")

file1<- file.choose()
file2<- file.choose()
file3<- file.choose()
file4 <- file.choose()
obj1 <- fasta2genlight(file1, chunk=6)
obj2 <- fasta2genlight(file2, chunk=6)
obj3 <- fasta2genlight(file3, chunk=6)
obj4 <- fasta2genlight(file4, chunk=6)

## look at extracted information
position(obj)
alleles(obj)
locNames(obj)

## plot positions of polymorphic sites
par(mfrow = c(2,2))
temp <- density(position(obj3), bw=10)
plot(temp, xlab="Position in the alignment", lwd=2, main="Location of the SNPs")
points(position(obj3), rep(0, nLoc(obj3)), pch="|", col="red")
