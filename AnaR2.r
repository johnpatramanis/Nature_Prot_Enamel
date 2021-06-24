
args = commandArgs(trailingOnly=TRUE)


library(ShortRead)

directory<-args[1]

setwd(directory)

genes<-readLines("Genes.txt")
sam<-args[2]

# 0 sample name
# 0.1 gene name
# 1. total number of sites in aln
# 2. total number of non-missing sites in the ancient sample
# 3. number of seg. sites in aln
# 4. number of seg. sites in the ancient sample
# 5. number of seg.&non-singletons sites in aln (here singleton means a seg site where only one sample has a difference)
# 6. number of seg-&non-singletons sites in the ancient sample (here singleton means a seg site where only one sample has a difference)
# 7. number of seg. sites where the ancient sample has a unique site



#################################################################
##SET UP FUNCTIONS

GetNotInfoPos<-function(x){
x<-x[x=="X" | x=="-"]
return(length(x))
}

GetSegSites<-function(x){
return(length(unique(x[as.character(x)!="-" & as.character(x)!="X"])))
}


GetSingletons<-function(x){
if(length(unique(x[as.character(x)!="-" & as.character(x)!="X"]))==1){
    return(0)
}else if(length(unique(x[as.character(x)!="-" & as.character(x)!="X"]))==2){
    if((table(x[as.character(x)!="-" & as.character(x)!="X"])[1]==1) |(table(x[as.character(x)!="-" & as.character(x)!="X"])[2]==1)){
        return(1)
    }else{
        return(0)
}
}else{
return(0)
}
}


GetAncientUnique<-function(x, aid){
    if(length(unique(x))>1){
        if(sum(x==x[aid])==1){
             return(1)
        }else{
              return(0)
        }
    }else{
          return(0)
    }
}


###############################################################################
###### FUNCTIONS END

###### GET INFO


tab<-NULL
for(g in 1:length(genes)){
        setwd(paste0(directory, genes[g]))
        f<-paste0(genes[g], "_aln.fa")
        
        fa<-readAAStringSet(f)
        print(fa)
        fatabble<-as.matrix(t(as.data.frame(strsplit(as.character(fa), ""))))
        print(fatabble)
        fanonmissing<-fatabble[,(fatabble[grep(sam, names(fa)),]!="-" & fatabble[grep(sam, names(fa)),]!="X")]
        print(fanonmissing)
        if(dim(fanonmissing)[2]>0){
            ## total sites
            TotalSites<-dim(fatabble)[2]
            SitesAncient<-dim(fanonmissing)[2]
            
            
            ##seg sites
            ss<-apply(fatabble, 2, GetSegSites)
            SegSites<-length(ss[ss>1])
            ss<-apply(fanonmissing, 2, GetSegSites)
            SegSitesAncient<-length(ss[ss>1])
            
            #singletons
            ss<-apply(fatabble, 2, GetSingletons)
            NonSingSites<-SegSites-length(ss[ss==1])
            ss<-apply(fanonmissing, 2, GetSingletons)
            NonSignAncient <-SegSitesAncient-length(ss[ss==1])
            
            #ancient unique
            AncientUnique<-sum(apply(fanonmissing, 2, GetAncientUnique, grep(sam, names(fa)))==1)
            tab<-rbind(tab, c(sam, genes[g], TotalSites, SitesAncient, SegSites, SegSitesAncient, NonSingSites, NonSignAncient, AncientUnique, f))
            
            
        }
        else{
            message("Error")
            }
    }

setwd(directory)

colnames(tab)<-c("Sample_name", "Gene_name", "Total_sites", "Sites_in_ancient", "Segregating_sites", "a_Segregating_sites", "Non_singletons", "a_Non_singletons", "Unique_ancient_sites", "File_name")

#change name here:
write.table(tab, file=paste0("Alignments_info_2021_","Sample_",sam,".txt"), col.names=T, row.names=F, quote=F, sep="\t")

#######################################################################################################################################################################################################################################


