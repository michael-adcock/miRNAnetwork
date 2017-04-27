args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Two arguments must be supplied. Input file and output file", call.=FALSE)
}

library("biclust", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
#setwd("~/Dropbox/project_1/data/NetworkMethod")

matrix_file <- read.csv(args[1], row.names=1)
matrix <- data.matrix(matrix_file, rownames.force = TRUE)

#if (args[3] == "T") {
#  matrix <- t(matrix)
#}

#print(dim(matrix))
clust1 <- biclust(x=matrix, method=BCrepBimax(), minr=2, minc=2, number=500)
print(clust1)
summary(clust1)
writeBiclusterResults(args[2], clust1, "ClusterMethod.py parameters: method=BCrepBimax(), minr=2, minc=2, number=500", rownames(matrix), colnames(matrix))
