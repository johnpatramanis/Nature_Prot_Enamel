###############################
# Rscript DataR3.r  
###############################

args<-(commandArgs(TRUE))

blastout<-"Neanderthal_ENAM_spliced.blast" 
geneid<-"ENAM"
sample_name<-"Neanderthal"


blastout<-args[1]
geneid<-args[2]
sample_name<-args[3]

library(ShortRead)

fout<-paste0(sample_name, "_", geneid, "_translated.fa")  #name of output fasta (protein)

con<-pipe(paste0("grep -e \"Identities\" -e \"Query\" -e \"Sbjct\" -e \"Length of query\" ", blastout, " | grep -v \"Query=\""))  #create pipe connection object for file
a<-scan(con, what="", sep="\n") #scan file
close(con) #
con<-pipe(paste0("grep \"letters)\" ", blastout)) #new pipe
len<-scan(con, what="", sep="\n")  #sequence
close(con)
len<-as.numeric(strsplit(strsplit(len, "\\(")[[1]][2], " ")[[1]][1]) #length of sequence
separator<-grep("Identities", a) #

if(length(separator)>1){
	a<-a[(separator[1]+1):(separator[2]-1)]
	b<-a[grep("Query",a)]
	tot<-as.numeric(strsplit(b[length(b)], " ")[[1]][length(strsplit(b[length(b)], " ")[[1]])])
	b<-as.numeric(strsplit(b[1], " ")[[1]][2])
	a<-a[grep("Sbjct",a)]
	a<-gsub("Sbjct: (\\d+)( *+)", "", a, perl=TRUE)
	a<-gsub(" (\\d+)", "", a, perl=TRUE)
	if(b!=1){
		a<-c(paste(collapse="", rep("X", b-1)), a)
	}
	if(tot<len){
		a<-c(a, paste(collapse="", rep("X", len-tot)))
	}
	a<-paste(collapse="", a)
	newseq<-AAStringSet(a)
	names(newseq)<-paste0(sample_name, "_", geneid)
	writeXStringSet(newseq, fout)
    
}else{ #if separator=0
	if(length(a)==1|length(a)==0){ #if no results at all! # seems to be called only for samples that lack AMELY, good!
		a=a[1] 
		print(paste0("No Results for 1 BLAST ", sample_name))            
	}else{
            a<-a[(separator[1]+1):(length(a))]
            b<-a[grep("Query",a,ignore.case=TRUE)]
            tot<-as.numeric(strsplit(b[length(b)], " ")[[1]][length(strsplit(b[length(b)], " ")[[1]])])
            b<-as.numeric(strsplit(b[1], " ")[[1]][2])
            a<-a[grep("Sbjct",a)]
            a<-gsub("Sbjct: (\\d+)( *+)", "", a, perl=TRUE)
            a<-gsub(" (\\d+)", "", a, perl=TRUE)
            if(b!=1){
                a<-c(paste(collapse="", rep("X", b-1)), a)
            }
            if(tot<len){
                a<-c(a, paste(collapse="", rep("X", len-tot)))
            }
            a<-paste(collapse="", a)
            newseq<-AAStringSet(a)
            names(newseq)<-paste0(sample_name, "_", geneid)
            writeXStringSet(newseq, fout)
	}        
}
