#!/bin/bash

################################################
###### Complete Bashcript, includes Rscripts
#######


################################################
####### REQUIREMENTS TO RUN:
#
# UNIX OS
#
# R (tested on version )
# angsd
# R packages: "ShortRead"
# BLAST tools (Blastn, Blastall
#


################################################
### Step 0
# Set up of file paths


MAIN_DIR='/path/to/main/directory/';

BAM_FILES='/path/to/where/bam_files_are_located/';

VCF_FILES='/path/to/where/vcf_files_are_located/';

GENE_LOCS='/path/to/Gene_locs.txt'; # text file with information on each genes location

STARTS='/path/to/starts.txt'; #text file with information on gene starts

PROT_REF_FILES='/path/to//Genes_human/'; #Folder with reference fasta files for each protein

EIT='/path/to/EIT/'; #Folder with Introns and Exons tables






################################################
### Step 1

#for bam files

cd $BAM_FILES;
for i in ./*bam;
do 
    pop=$(basename $i |cut -f 1 -d "_"); 
    sample=$(basename $i |cut -f 2 -d "_");    
    chr=$(basename $i  |cut -d "_" -f 3 |perl -pe 's/.bam//g;');    
    angsd -minQ 20 -minMapQ 30 -doFasta 2 -doCounts 1 -basesPerLine 60 -i $i -r chr$chr -out $pop"_"$sample"_"$chr"";
done







################################################
### Step 2

>#for VCF files
>cd $VCF_FILES;
>for i in *vcf.gz;
>    do
>        chr=$(echo $i |cut -d "_" -f 1);
>      na=$(echo $i |cut -d "." -f 1);
>     for GRP in "Pongo" "Gorilla_beringei" "Gorilla_gorilla" "Homo_sapiens" "Pan_paniscus" "Pan_troglodytes-" "Pan_troglodytes_ellioti" "Pan_troglodytes_schweinfurthii" "Pan_troglodytes_troglodytes" "Pan_troglodytes_verus";
>            do
>               cat SampleNames.txt | grep $GRP |while read line;
>                do
>                    sample=$(echo $line |cut -d " " -f 1);
>                    ids=$(echo $line |cut -d " " -f 2);
>                  echo $na $chr $sample $ids;
>                  zcat $i | cut -f 1-9,"$ids" > "$sample"_"$na"_TMP.vcf ;
>                  perl ./CnsFromVCF_prots.pl -i "$sample"_"$na"_TMP.vcf -s "$sample" -chr "$chr";
>                  rm "$sample"_"$na"_TMP.vcf;
>              done;
>           done;
>    done;





################################################
### Step 3
# Rfiles should be in the same Dir

>cd $BAM_FILES;
>Rscript DataR1.r $GENE_LOCS;
>cd $VCF_FILES;
>Rscript DataR1.r $GENE_LOCS;









################################################
### Step 4
# Rfiles should be in the same Dir

>cd $BAM_FILES;
>Rscript DataR2.r $STARTS $EIT
>cd $VCF_FILES;
>Rscript DataR2.r $STARTS $EIT





################################################
### Step 5
#for Bam files

>cd $BAM_FILES

>for i in *_spliced.fa; do makeblastdb -dbtype nucl -in $i; done ##creates a BLAST database for each file!

>ls *spliced.fa |cut -f 1,2 -d "_" |sort |uniq > SAMPLES 

>cat $GENE_LOCS |while read line;
>    do gene=$(echo $line |cut -d " " -f 1);
>       fref=$(echo $line |cut -d " " -f 5);
>      cat SAMPLES |while read sams;
> 	 do blastall  -p tblastn -i $PROT_REF_FILES"$fref"  -d "$sams"_"$gene"_spliced.fa -o "$sams"_"$gene"_spliced.blast -F F -E 32767 -G 32767 -n T -m 0 -M PAM70;
>          done;
>    done;

>Rscript DataR3.r 


#for VCF files files


>cd $VCF_FILES;

>for i in *_spliced.fa; do makeblastdb -dbtype nucl -in $i; done ##creates a BLAST database for each file!

>ls *spliced.fa |cut -f 1,2 -d "_" |sort |uniq > SAMPLES 

>cat $GENE_LOCS |while read line;
>    do gene=$(echo $line |cut -d " " -f 1);
>       fref=$(echo $line |cut -d " " -f 5);
>      cat SAMPLES |while read sams;
> 	 do blastall  -p tblastn -i $PROT_REF_FILES"$fref"  -d "$sams"_"$gene"_spliced.fa -o "$sams"_"$gene"_spliced.blast -F F -E 32767 -G 32767 -n T -m 0 -M PAM70;
>          done;
>    done;

>Rscript DataR3.r 




################################################
### Step 6

>cd $BAM_FILES
>Rscript DataR4.r
>cd $VCF_FILES
>Rscript DataR4.r







################################################
### Step 7

>cd $MAIN_DIR

>mkdir REFERENCE_DATASET #create dir

>rm REFERENCE_DATASET/*_PROT_REFERENCE.fa   #empty it, if it exists

>cd $BAM_FILES
>ls *_translated.fa |cut -f 3 -d "_" |sort |uniq |while read line; do touch $MAIN_DIR/REFERENCE_DATASET/"$line"_PROT_REFERENCE.fa; cat *"$line"_translated.fa >$MAIN_DIR/REFERENCE_DATASET/"$line"_PROT_REFERENCE.fa; done;  

>touch $MAIN_DIR/REFERENCE_DATASET/ALL_PROT_REFERENCE.fa
>cat *_translated.fa >ALL_PROT_REFERENCE.fa;
>mv ALL_PROT_REFERENCE.fa ../REFERENCE_DATASET/ALL_PROT_REFERENCE.fa;

#For vcf as well...

>cd $VCF_FILES
>ls *_translated.fa |cut -f 3 -d "_" |sort |uniq |while read line; do cat *"$line"_translated.fa >>$MAIN_DIR/REFERENCE_DATASET/"$line"_PROT_REFERENCE.fa; done
>cat *_translated.fa >>../REFERENCE_DATASET/ALL_PROT_REFERENCE.fa;




