# Pathway Enrichment Test
Conduct functional analysis using the hypergeometric test to find any enrichment among pathways.

Example inputs and outputs have been provided in the corresponding folders.

## Two inputs will be prompted
### 1. The gene pathway file
A tab deliminated file that contains two columns. The first column contains the gene id and the second column contains the pathway the gene is in. Each relationship should be seperated by a newline, so a gene may appear on multiple lines to show its relationsihp to each pathway.
### 2. A gene of interest file
A one-column file with each gene id seperated by a newline. The gene ids used in this file must correspond with the gene ids provided in the gene pathway file. An example is a list of differentially expressed genes.

## Output
Table with the following columns will be outputted:
1. pathway - name of the pathway
2. wht.drawn - number of gene of interest in pathway
3. wht.in.urn - number of genes in pathway
4. black.in.urn - total genes mapped to other pathways
5. total.draw - total number of genes of interest that belong to a pathway
6. p.value - p-value as a result of the hypergeometric test
7. p.adjust - adjusted p-value using the `fdr` method within R
