library("readxl")
library(readr)
library(GSA)
library("xlsx")


filename1 <- 'gprofiler_full_mmusculus.name.gmt'
ugs <- GSA.read.gmt(filename1) 
universal_gene_set=unlist(ugs$genesets)   ## Vector


df <- data.frame(matrix(unlist(universal_gene_set), nrow=length(universal_gene_set), byrow=TRUE))
write_delim(df, "Universal_genes_Final_from_gprofiler.txt", delim = "\n", append=TRUE)





