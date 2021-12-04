library(biomaRt)
library(readr)
ensembl <- useMart("ensembl")
ensembl_mart<-useMart("ensembl",dataset= "mmusculus_gene_ensembl")  

#mouse: mmusculus_gene_ensembl

#human: hsapiens_gene_ensembl 


for(i in 1:19)
{
  #chromosome <- getBM(attributes = ("hgnc_symbol"), filter="chromosome_name", values=i, mart= ensembl_mart) #Human
  #chromosome <- getBM(attributes = ("wikigene_name"), filter="chromosome_name", values=i, mart= ensembl_mart)
  chromosome <- getBM(attributes = ("external_gene_name"), filter="chromosome_name", values=i, mart= ensembl_mart)
  write_delim(chromosome, "All_genes_Final.txt", delim = "\n", append=TRUE)
}

#chromosome <- getBM(attributes = ("hgnc_symbol"), filter="chromosome_name", values="x", mart= ensembl_mart) #Human
chromosome <- getBM(attributes = ("external_gene_name"), filter="chromosome_name", values="x", mart= ensembl_mart)
write_delim(chromosome, "All_genes_Final.txt", delim = "\n", append=TRUE)

#chromosome <- getBM(attributes = ("hgnc_symbol"), filter="chromosome_name", values="y", mart= ensembl_mart) #Human
chromosome <- getBM(attributes = ("external_gene_name"), filter="chromosome_name", values="y", mart= ensembl_mart)
write_delim(chromosome, "All_genes_Final.txt", delim = "\n", append=TRUE)



#chromosome <- getBM(attributes = ("wikigene_name"), mart= ensembl_mart)
#write_delim(chromosome, "All_genes_Final.txt", delim = "\n", append=TRUE)