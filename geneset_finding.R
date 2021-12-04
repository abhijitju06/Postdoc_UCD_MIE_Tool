library("readxl")
library(readr)
my_data <- read_excel("MousePath_All.xlsx", sheet = "MousePath")
s=length(my_data$genesSym)
for (i in 1:s)
{
  genes <- my_data$genesSym[i]
  genes=gsub(",","\n",genes)
  df <- data.frame(matrix(unlist(genes), nrow=length(genes), byrow=TRUE))
  write_delim(df, "New_All_genes_Final.txt", delim = "\n", append=TRUE)
}

