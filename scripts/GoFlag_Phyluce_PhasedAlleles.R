library(ape)
library(phytools)
appTrees <- read.tree(file.choose())

View(appTrees)
appMatrix <- matrix(data = NA, ncol = 8, nrow = length(appTrees)*2)
colnames(appMatrix) <- c("eng-app0", "val-app0", "app0-app1", "val-app1", "eng-app1", "app0-app1", "App0 Min Distance", "App1 Min Distance")
j <- 1
for (i in seq(1,length(appTrees),1)) {
  appMatrix[j,1] <- fastDist(appTrees[[i]], "Isoetes_engelmannii", "Isoetes_appPhase0")
  appMatrix[j,2] <- fastDist(appTrees[[i]], "Isoetes_valida", "Isoetes_appPhase0")
  appMatrix[j,3] <- fastDist(appTrees[[i]], "Isoetes_appPhase0", "Isoetes_appPhase1")
  appMatrix[j+1,4] <- fastDist(appTrees[[i]], "Isoetes_valida", "Isoetes_appPhase1")
  appMatrix[j+1,5] <- fastDist(appTrees[[i]], "Isoetes_engelmannii", "Isoetes_appPhase1")
  appMatrix[j+1,6] <- fastDist(appTrees[[i]], "Isoetes_appPhase0", "Isoetes_appPhase1")
  appMatrix[j,7] <- colnames(appMatrix)[which.min(appMatrix[i,])]
  appMatrix[j+1,8] <- colnames(appMatrix)[which.min(appMatrix[i+1,])]
  j <- j + 2
}

View(appMatrix)

phaseSums <- matrix(data = NA, nrow = 1, ncol = 5)
colnames(phaseSums) <- c("eng-app0", "val-app0", "eng-app1", "val-app1", "app0-app1")
phaseSums[1,1] <- sum(appMatrix[,7] == "eng-app0", na.rm = TRUE)
phaseSums[1,2] <- sum(appMatrix[,7] == "val-app0", na.rm = TRUE)
phaseSums[1,3] <- sum(appMatrix[,8] == "eng-app1", na.rm = TRUE)
phaseSums[1,4] <- sum(appMatrix[,8] == "val-app1", na.rm = TRUE)
phaseSums[1,5] <- sum(appMatrix[,8] == "app0-app1", na.rm = TRUE)

View(phaseSums)
x <- barplot(phaseSums, xaxt = "n", ylab = "Number of genes supporting relationship")
labs <- colnames(phaseSums)
text(cex=1, x=x-.25, y=-1, labs, xpd=TRUE, srt=45)
