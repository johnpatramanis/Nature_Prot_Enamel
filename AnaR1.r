args = commandArgs(trailingOnly=TRUE)




### If you need to install ShortRead, check the instructions here: https://bioconductor.org/packages/release/bioc/html/ShortRead.html
library(ShortRead)

directory=args[1]

# This correspond to the protein sequences for the reference samples
fa<-readAAStringSet(args[2])


genes<-readLines("Genes.txt")



for(i in 1:length(genes)){
d<-paste0(directory, genes[i])
setwd(d)
curfa<-fa[grep(genes[i], names(fa))]  ### biostring with only one gene (all fasta seqs that have the gene name in their name), all sample sequences tied with their sample name


writeXStringSet(curfa, paste0(genes[i], "_o.fa"))
}
