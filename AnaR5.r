
args = commandArgs(trailingOnly=TRUE)



library("phyclust")

directory<-args[1]
setwd(directory)

genes<-readLines("Genes.txt")

for(g in 1:length(genes)){
    setwd(paste0(directory,genes[g]))
    fasta_data=read.fasta(paste0("./",genes[g],"_aln_e.fa"), code.type ="AMINO_ACID",aligned=TRUE)
    names=fasta_data$seqname
    for (i in 1:length(names)){
        names[i]=paste0(substring(names[i], 1, nchar(names[i])),"\t")
    }


    write.phylip( fasta_data$org , paste0("./",genes[g],"_aln_e.phy") ,code.type = "AMINO_ACID",seqname=names, width.seqname =100)

# backup method if errors
# write.phylip.format( fasta_data$org.code , paste0("./",genes[g],"_aln_e.phy"),seqname=names, width.seqname =10)
}

######### Turn this one into Phyl as well


setwd(directory)
concat_fasta=read.fasta("CONCATINATED_o.fa", code.type ="AMINO_ACID",aligned=TRUE)

names=concat_fasta$seqname 
for (i in 1:length(names)){
    names[i]=paste0(substring(names[i], 1, nchar(names[i])),"\t")
}

write.phylip( concat_fasta$org , "./CONCATINATED_aln_e.phy" ,code.type = "AMINO_ACID",seqname=names, width.seqname =100)
