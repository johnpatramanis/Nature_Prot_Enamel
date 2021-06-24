print("Moving to I,L switching")
#######################################################################################################################################################################################################################################
#####Switch Animo Acids I and L

library(ShortRead)

setwd(directory)



genes<-readLines("Genes.txt")

#### should be changed per sample 

for(g in 1:length(genes)){              #### for every gene
    setwd(paste0(directory, genes[g]))   ### go to the dir
    f<-paste0(genes[g], "_aln.fa")       ## grab the aligned fasta file
    fa<-readAAStringSet(f)               ## get it as a biostring
    fatabble<-as.matrix(t(as.data.frame(strsplit(as.character(fa), ""))))  # transform it to a matrix
    aid<-grep(sam, names(fa))
    Jsites<-which(fatabble[aid,]=="I" | fatabble[aid,]=="L")    # find which sites of sample have either an I or L
    if(length(Jsites)>0){   #for every one of these sites
        for(s in 1:length(Jsites)){
            cursite<-fatabble[,Jsites[s]]
            optssite<-as.character(cursite[-grep(sam, names(cursite))])
            if(length(table(optssite[optssite=="L" | optssite=="I"]))==1){
                cursite[cursite=="L" | cursite=="I"]<-names(table(optssite[optssite=="L" | optssite=="I"]))
            }
            else{
                cursite[cursite=="L" | cursite=="I"]<-"L"
            }
            fatabble[,Jsites[s]]<-cursite
            }
        newseq<-apply(fatabble, 1, paste, collapse="")
        names(newseq)<-names(fa)
        writeXStringSet(AAStringSet(newseq), gsub(".fa", "_e.fa", f))
        }
    else{
        writeXStringSet(fa, gsub(".fa", "_e.fa", f))
    }
}