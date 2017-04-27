library(ReactomePA)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Three arguments must be supplied. Input file, output file and pvalue", call.=FALSE)
}

gene_ids_file <- file(args[1], open="r")
genes <-readLines(gene_ids_file)

output <- args[2]
if(file.exists(output)){
  file.remove(output)
}

pvalue <- as.numeric(args[3])
x <- enrichPathway(gene=genes,pvalueCutoff=pvalue, readable=T)
print(as.data.frame(x))
write.table(x, output, sep="\t")
