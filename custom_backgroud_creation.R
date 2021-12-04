###### Load Library ######
library(GSA)
######  MIE1  #############
filename1 <- 'MIE1.gmt'
gs1 <- GSA.read.gmt(filename1) ##List of Lists with gene-list as first list
MIE1_gene_set=unique(unlist(gs1[1]))   ## Vector
#####################################################################

######  MIE2  #############
filename2 <- 'MIE2.gmt'
gs2 <- GSA.read.gmt(filename2) ##List of Lists with gene-list as first list
MIE2_gene_set=unique(unlist(gs2[1]))   ## Vector
#####################################################################

######  MIE3  #############
filename3 <- 'MIE3.gmt'
gs3 <- GSA.read.gmt(filename3)  ##List of Lists with gene-list as first list
MIE3_gene_set=unique(unlist(gs3[1]))    ## Vector
#####################################################################

######A vector with all genes, including the genes of interest
universe_genes=union(MIE1_gene_set,MIE2_gene_set)
universe_genes=union(universe_genes,MIE3_gene_set)

write.table(universe_genes, file = "custom_background.txt", sep = "\n",
            row.names = FALSE, col.names = FALSE)