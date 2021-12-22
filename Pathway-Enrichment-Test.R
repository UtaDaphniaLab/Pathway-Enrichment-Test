library(dplyr)
library(tidyr)


# load gene pathway table
gene.pway.table = read.table(choose.files(caption="Select gene pathway file"), 
                             sep="\t", header=FALSE, col.names=c("gene", "pathway")) # choose gene-pathway-table
num.genes=nrow(gene.pway.table %>% count(gene)) # number of unique genes mapped to pathways
mapped.genes = length(gene.pway.table$gene) # instances of genes mapped to pathways (including genes repeated many times)

# calculates number of genes mapped to each pathway
pway.table = as.data.frame(gene.pway.table %>% count(pathway))
colnames(pway.table) = c("pathway", "total_genes")
npathway = length(pway.table$pathway) # unique pathways found in reference


# load gene list
goi.list = read.table(choose.files(caption="Select gene of interest interest file"), 
                      sep="\t", header=F, col.names=c("gene")) # choose gene of interest (goi) list
goi.mapped = merge (goi.list, gene.pway.table, by.x="gene", by.y="gene", all.x=FALSE, all.y=FALSE) # find pathway for each gene. remove genes not mapped to pathway
goi.pway.count = goi.mapped %>% count(pathway) # counts number of genes of interest in each pathway

# collects data for and runs hypergeometric test
hypergeometric.table = merge(goi.pway.count, pway.table, by.x= "pathway", by.y="pathway") # create table to hold data for hypergeometric test
colnames(hypergeometric.table) = c("pathway", "wht.drawn", "wht.in.urn")
# black in urn: number of genes not mapped to pathway of interest
hypergeometric.table$blk.in.urn = num.genes - hypergeometric.table$wht.in.urn
# finds the number of unique DEGs/genes of interest that are mapped to a pathway (total.draw)
hypergeometric.table$total.draw = nrow(goi.mapped %>% count(gene))
#perform hypergeometric test
p.value = apply(hypergeometric.table, 1, function(x) phyper(as.numeric(x[2])-1, as.numeric(x[3]), as.numeric(x[4]), as.numeric(x[5]),lower.tail= FALSE))
hypergeometric.table = data.frame(hypergeometric.table, p.value)

hypergeometric.table$p.adjust = round(p.adjust(hypergeometric.table$p.value, "fdr"),4)

# sorts by ascending p value
hypergeometric.table = hypergeometric.table[order(hypergeometric.table$p.value), ]
write.table(hypergeometric.table, file=file.choose(), sep = "\t", col.names=TRUE, row.names = FALSE, quote=FALSE)
