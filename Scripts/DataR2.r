###############################
# Rscript DataR2.r  
###############################

args<-(commandArgs(TRUE))

s<-"Data/starts.txt" 
exonb<-"Data/ENAM_ei.txt"
geneid<-"ENAM"
sample_name<-"Neanderthal"


s<-args[1]
exonb <-args[2]
geneid<-args[3]
sample_name<-args[4]

library(ShortRead)

# name of each gene/ where each gene starts / which strand
info<-read.table(s, h=T, as.is=T) 

fa<-readFasta(paste0(sample_name, "_", geneid, ".fa")) #readFasta loads the fasta file into R
fa<-DNAString(as.character(sread(fa))) # Turn Fasta into a string

exon_boundaries<-read.table(exonb, as.is=T, sep="\t") # Exon / Intron file
    
if(info[info[,1]==geneid,3]=="+"){ # If the gene of the fasta file is on the (+) strand
	starts<-info[info[,1]==geneid,2] #grab first start position
	ends<-NULL
	for(i in 1:length(exon_boundaries[,1])){   #use intron position to create chunks of introns
		ends<-c(ends, ((starts[length(starts)]+exon_boundaries[i,2])-1)) # get end position from last start position + length of intron/exon -1
		starts<-c(starts, (starts[length(starts)]+exon_boundaries[i,2])) # next start position is last end position (+1)
	}
	starts<-starts[1:length(exon_boundaries[,1])] # remove the last start position
        
	starts<-starts[-grep("Intron", exon_boundaries[,1])] # remove intron starts from list
	ends<-ends[-grep("Intron", exon_boundaries[,1])]    #remove intron ends from list

}else{ # if gene of the fasta file is on (-) strand, same as above but the reverse logic (move from right to left)
	ends<-info[info[,1]==geneid,2] # what previously would be start is here the end
	starts<-NULL #the same logic as above
	for(i in 1:length(exon_boundaries[,1])){ #the same logic as in the above loop, but we are moving to the left, so starts are bigger numbers than their ends
		starts<-c(starts, ends[length(ends)]-exon_boundaries[i,2]+1)
		ends<-c(ends, ends[length(ends)]-exon_boundaries[i,2])
	}
        
	ends<-ends[1:length(exon_boundaries[,1])]   #again remove last part
        #remove introns
	starts<-starts[-grep("Intron", exon_boundaries[,1])]
	ends<-ends[-grep("Intron", exon_boundaries[,1])]
	starts<-rev(starts) #flip them into canonical order, from smaller number to larger
	ends <-rev(ends)
}
    
seqs<-NULL
#now we use those starts / ends pairs to isolate the exons only from the sequence
for(i in 1:length(starts)){
	seqs<-paste(sep="", seqs, as.character(fa[starts[i]:ends[i]]))
}
seqs<-gsub(" ", "", seqs)
    
ids<-paste0(sample_name, "_", geneid, "_spliced")
    
newseq<-ShortRead(sread=DNAStringSet(seqs), id=BStringSet(ids))
    
writeFasta(newseq, paste0(sample_name, "_", geneid, "_spliced.fa"))

