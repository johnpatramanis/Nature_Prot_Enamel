args = commandArgs(trailingOnly=TRUE)
library(ShortRead)

f<-dir(pattern="\\.fa.gz$")   #thsoe are created from previous step, or from other pipelines
genes<-read.table(args[1], as.is=T) # this file contains the chromosome/positions for the genes of interest

for(s in 1:length(f)){
    fa<-readFasta(f[s])
    chr<-gsub(".fa.gz", "", sapply(strsplit(basename(f[s]), "_"), "[[", 3))
	chr<-
    sam<-gsub(".fa.gz", "", sapply(strsplit(basename(f[s]), "_"), "[[", 2))
    pop<-gsub(".fa.gz", "", sapply(strsplit(basename(f[s]), "_"), "[[", 1))
    curgenes<-genes[genes[,2]==chr,]

    for(i in 1:length(curgenes[,1])){
        seq<-as.character(DNAString(as.character(sread(fa)))[curgenes[i,3]:curgenes[i,4]]) ## Isolate location of gene in chromosome, start-stop
        newseq<-ShortRead(sread=DNAStringSet(seq), id=BStringSet(paste0(pop, "_", sam, "_", curgenes[i,1], "_", curgenes[i,2], "_", curgenes[i,3], "_", curgenes[i,4]))) # prepare fasta sequence
        writeFasta(newseq, paste0(pop, "_", sam, "_", curgenes[i,1], ".fa")) # write it out
    }
}

