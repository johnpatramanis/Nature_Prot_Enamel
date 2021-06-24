library(ShortRead)

f<-dir(pattern="\\_spliced.blast")    #load all blast files
genes<-gsub("_spliced.blast", "", f)     #get the names for each blast file

fout<-paste0(genes, "_translated.fa")  #name of output fasta (protein)
h<-genes
h<-paste(sep="_", sapply(strsplit(h, "_"), "[[", 1), sapply(strsplit(h, "_"), "[[", 2), sapply(strsplit(h, "_"), "[[", 3)) # just the sample name

for(i in 1:length(f)){ #for each blast file
    blout<-f[i]   #grab each file
    zz<-pipe(paste0("grep -e \"Identities\" -e \"Query\" -e \"Sbjct\" -e \"Length of query\" ", blout, " | grep -v \"Query=\""))  #create pipe connection object for file
    a<-scan(zz, what="", sep="\n") #scan pipe ?
    close(zz) #
    zz<-pipe(paste0("grep \"letters)\" ", blout)) #new pipe
    len<-scan(zz, what="", sep="\n")  #sequence
    close(zz)
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
        names(newseq)<-h[i]
        writeXStringSet(newseq, fout[i])
    
    }else{#if separator=0
        if(length(a)==1|length(a)==0){ #if no results at all! # seems to be called only for samples that lack AMELY, good!
            a=a[1] 
            print(paste0("No Results for 1 BLAST ",h[i]))
            # b<-a[grep("Query",a,ignore.case=TRUE)]
            # tot<-as.numeric(strsplit(b[length(b)], " ")[[1]][length(strsplit(b[length(b)], " ")[[1]])])
            # a=rep("X", tot)
            # a<-paste(collapse="", a)            
            # newseq<-AAStringSet(a)
            # names(newseq)<-h[i]
            # writeXStringSet(newseq, fout[i])
            
            
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
            names(newseq)<-h[i]
            writeXStringSet(newseq, fout[i])
            }        
    }
}
