###############################
# Rscript AnaR2.r  
#
# Switch Animo Acids I and L
###############################

args<-(commandArgs(TRUE))

alignment<-"ENAM_aligned.fa" 
geneid<-"ENAM"
sample_name<-"H_antecessor"

alingment<-args[1]
geneid<-args[2]
sample_name<-args[3]

library(ShortRead)

#

fa<-readAAStringSet(alignment)
fatabble<-as.matrix(t(as.data.frame(strsplit(as.character(fa), ""))))
aid<-grep(sample_name, names(fa))
Jsites<-which(fatabble[aid,]=="I" | fatabble[aid,]=="L")    # find which sites of sample have either an I or L
if(length(Jsites)>0){   #for every one of these sites
	count_cs<-0
	for(s in 1:length(Jsites)){
		cursite<-fatabble[,Jsites[s]]
		optssite<-as.character(cursite[-grep(sample_name, names(cursite))])
		if(length(table(optssite[optssite=="L" | optssite=="I"]))==1){
			cursite[cursite=="L" | cursite=="I"]<-names(table(optssite[optssite=="L" | optssite=="I"]))
			count_cs<-count_cs+1
		}else{
			# decide what to do on this:
			cursite[cursite=="L" | cursite=="I"]<-"L"
			#cursite[grep(sample_name, names(cursite))]<-"X"
			count_cs<-count_cs+1
		}
		fatabble[,Jsites[s]]<-cursite
	}
	newseq<-apply(fatabble, 1, paste, collapse="")
	names(newseq)<-names(fa)
	writeXStringSet(AAStringSet(newseq), gsub(".fa", "_noIso.fa", alignment))
}else{
        writeXStringSet(fa, gsub(".fa", "_noIso.fa", alignment))
}
message(paste0("There were ", count_cs, " sites changed."))

