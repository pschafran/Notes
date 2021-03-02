csv <- read.table("~/Desktop/binmatrix.csv", sep=",")
binmatrix <- matrix(as.numeric(csv), ncol=100)
image(binmatrix, col= collist, xaxt="n", yaxt="n")