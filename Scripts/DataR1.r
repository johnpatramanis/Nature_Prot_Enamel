###############################
# Rscript DataR1.r  
###############################

args<-(commandArgs(TRUE))

#genepos<-"Data/Gene_locs.txt" 
#fafile<-"Neanderthal_ENAM.fa.gz"
#geneid<-"ENAM"
#sample_name<-"Neanderthal"

genepos<-args[1]
fafile<-args[2]
geneid<-args[3]
sample_name<-args[4]

library(ShortRead)

genes<-read.table(genepos, as.is=T)

fa<-readFasta(fafile)

# Isolate location of gene in chromosome, start-stop
seq<-as.character(DNAString(as.character(sread(fa)))[genes[genes[,1]==geneid,3]:genes[genes[,1]==geneid,4]]) 

# prepare fasta sequence
newseq<-ShortRead(sread=DNAStringSet(seq), id=BStringSet(paste0(sample_name, "_", geneid))) 

writeFasta(newseq, paste0(sample_name, "_", geneid, ".fa"))

