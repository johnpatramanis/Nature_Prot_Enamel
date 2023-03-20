###############################
# Rscript AnaR3.r  
#
# Switch Animo Acids I and L
###############################

args<-(commandArgs(TRUE))

alignment<-"ENAM_aligned_noIso.fa" 
outfile<-"ENAM_aligned_noIso.phy"

alingment<-args[1]
outfile<-args[2]

library("phyclust")

fa<-read.fasta(alignment)

mat<-fa$org.code
rownames(mat)<-fa$seqname

write.phylip.format(mat, filename=outfile)

#