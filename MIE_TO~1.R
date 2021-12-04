#' Over representation test via hyper-geometric distribution
#'
#' @param signature A vector of symbols
#' @param genesets A list of gene sets
#' @param background Size of background population genes
#' @return A list of data 
#'
#' @importFrom stats phyper p.adjust
#' @keywords internal
hyper_enrichment <- function(signature,
                              genesets,custom_background, genesets_description,
                              background=length(custom_background)) {

    if (!is(signature, "vector")) stop("Expected signature to be a vector of symbols\n")
    if (!is(genesets, "list")) stop("Expected genesets to be a list of genesets\n")
    
    signature <- unique(signature)
    genesets <- lapply(genesets, unique)
    
    # Construct table
    signature.found <- signature[signature %in% unique(unlist(genesets))]
    n.hits <- sapply(genesets, function(x, y) length(intersect(x, y)), signature.found)
    n.drawn <- length(signature)
    n.genesets <- sapply(genesets, length)
    n.left <- background-n.genesets
    
    # Hypergeometric test
        pvals <- suppressWarnings(stats::phyper(q=n.hits-1,
                                                m=n.genesets,
                                                n=n.left,
                                                k=n.drawn,
                                                lower.tail=FALSE))
        suppressWarnings(if(is.nan(pvals))pvals=1)
    
    
    
    # Format data
    data <- data.frame(label=names(genesets),
                       label=names(genesets_description),
                       pval=signif(stats::p.adjust(pvals, method="fdr"), 3),
                       signature=length(signature),
                       geneset=n.genesets,
                       overlap=n.hits,
                       background=background,
                       hits=sapply(genesets, function(x, y) paste(intersect(x, y), collapse=','), signature.found),
                       stringsAsFactors=FALSE)
    
    
    return(list(data=data))
}


############# CALCULATION #########################################################################
library("xlsx")
library(GSA)
###################################################################################################
######  MIE1  #############
filename1 <- 'MIE1.gmt'
gs1 <- GSA.read.gmt(filename1) ##List of Lists with gene-list as first list
k1 <- length(gs1$genesets)
genesets1 <- list()
genesets1_description <-list()
for (i in 1:k1)
{
    temp <- gs1$genesets[i]
    names(temp)<-gs1$geneset.names[i]
    genesets1<- c(genesets1,temp)
    names(temp) <- gs1$geneset.descriptions[i]
    genesets1_description<- c(genesets1_description,temp)
}

######  MIE2  #############
filename2 <- 'MIE2.gmt'
gs2 <- GSA.read.gmt(filename2) ##List of Lists with gene-list as first list
k2 <- length(gs2$genesets)
genesets2 <- list()
genesets2_description <-list()
for (i in 1:k2)
{
    temp <- gs2$genesets[i]
    names(temp)<-gs2$geneset.names[i]
    genesets2<- c(genesets2,temp)
    names(temp) <- gs2$geneset.descriptions[i]
    genesets2_description<- c(genesets2_description,temp)
}
######  MIE3  #############
filename3 <- 'MIE3.gmt'
gs3 <- GSA.read.gmt(filename3)  ##List of Lists with gene-list as first list
k3 <- length(gs3$genesets)
genesets3 <- list()
genesets3_description <-list()
for (i in 1:k3)
{
    temp <- gs3$genesets[i]
    names(temp)<-gs3$geneset.names[i]
    genesets3<- c(genesets3,temp)
    names(temp) <- gs3$geneset.descriptions[i]
    genesets3_description<- c(genesets3_description,temp)
}

##########################################################################################################
######## Background ######################################################################################
custom_background <-read.delim('All_genes_Final.txt',header=FALSE,sep ="\n")
custom_background <-unlist(custom_background) ##Vector

signature <-read.delim('DEGlistsimple.txt',header=FALSE,sep ="\n") ##List
signature=unique(unlist(signature)) ##Vector
signature=intersect(signature,custom_background)

##########################################################################################################

wb<-createWorkbook(type="xlsx")
heading<-list("Geneset_name","Event","P-Value","No_Signature","No_Genes","overlap","Background","Hits")
###########################################################################################################
################################# MIE1 ####################################################################
ws1 <- createSheet(wb, sheetName = "MIE1")
setColumnWidth(ws1, colIndex=1:1, colWidth= 25)
setColumnWidth(ws1, colIndex=2:2, colWidth= 25)
setColumnWidth(ws1, colIndex=3:3, colWidth= 25)
setColumnWidth(ws1, colIndex=4:4, colWidth= 25)
setColumnWidth(ws1, colIndex=5:5, colWidth= 25)
setColumnWidth(ws1, colIndex=6:6, colWidth= 25)
setColumnWidth(ws1, colIndex=7:7, colWidth= 25)
setColumnWidth(ws1, colIndex=8:8, colWidth= 25)
setColumnWidth(ws1, colIndex=9:9, colWidth= 255)
addDataFrame(heading,sheet = ws1,row.names=FALSE, col.names=FALSE,startRow=1)

startRow = 2
result1 <- hyper_enrichment(signature, genesets1,custom_background, genesets1_description)

n=length(result1$data$label)
for (i in 1:n)
{
    if(result1$data[i,]$pval<0.05)
    {
        addDataFrame(result1$data[i,], sheet = ws1,row.names=FALSE, col.names=FALSE,startRow=startRow)
        startRow = startRow + nrow(result1$data[i,])
    }
    
}


###########################################################################################################
################################# MIE2 ####################################################################
ws2 <- createSheet(wb, sheetName = "MIE2")
setColumnWidth(ws2, colIndex=1:1, colWidth= 25)
setColumnWidth(ws2, colIndex=2:2, colWidth= 25)
setColumnWidth(ws2, colIndex=3:3, colWidth= 25)
setColumnWidth(ws2, colIndex=4:4, colWidth= 25)
setColumnWidth(ws2, colIndex=5:5, colWidth= 25)
setColumnWidth(ws2, colIndex=6:6, colWidth= 25)
setColumnWidth(ws2, colIndex=7:7, colWidth= 25)
setColumnWidth(ws2, colIndex=8:8, colWidth= 255)
addDataFrame(heading,sheet = ws2,row.names=FALSE, col.names=FALSE,startRow=1)

startRow = 2
result2 <- hyper_enrichment(signature, genesets2,custom_background, genesets2_description)

n=length(result2$data$label)
for (i in 1:n)
{
    if(result2$data[i,]$pval<0.05)
    {
        addDataFrame(result2$data[i,], sheet = ws2,row.names=FALSE, col.names=FALSE,startRow=startRow)
        startRow = startRow + nrow(result2$data[i,])
    }
    
}

###########################################################################################################
################################# MIE3 ####################################################################
ws3 <- createSheet(wb, sheetName = "MIE3")
setColumnWidth(ws3, colIndex=1:1, colWidth= 25)
setColumnWidth(ws3, colIndex=2:2, colWidth= 25)
setColumnWidth(ws3, colIndex=3:3, colWidth= 25)
setColumnWidth(ws3, colIndex=4:4, colWidth= 25)
setColumnWidth(ws3, colIndex=5:5, colWidth= 25)
setColumnWidth(ws3, colIndex=6:6, colWidth= 25)
setColumnWidth(ws3, colIndex=7:7, colWidth= 25)
setColumnWidth(ws3, colIndex=8:8, colWidth= 255)
addDataFrame(heading,sheet = ws3,row.names=FALSE, col.names=FALSE,startRow=1)

startRow = 2
result3 <- hyper_enrichment(signature, genesets3,custom_background, genesets3_description)

n=length(result3$data$label)
for (i in 1:n)
{
    if(result3$data[i,]$pval<0.05)
    {
        addDataFrame(result3$data[i,], sheet = ws3,row.names=FALSE, col.names=FALSE,startRow=startRow)
        startRow = startRow + nrow(result3$data[i,])
    }
    
}

saveWorkbook(wb, file = "My_results_DEGlistsimple.xlsx")







