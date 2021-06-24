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
#
# R packages: "ShortRead", "stringr" , "data.table", "phyclust"
#
# mafft 
# OR CLUSTAL OMEGA
#
# PhyML

################################################
### Step 1
# Set up of file paths


HOME_DIR='/home/path/to/desired/directory/';
REF='REFERENCE_DATASET.fa'; #fasta file of reference dataset
SAMPLE='SAMPLE.fa'; #fasta file of all sample’s proteins
SAMPLE_NAME='SAMPLE'; #The name of the sample
COMPLETENESS_CUTOFF='0.10'; #Proteins under this 








################################################
### Step 2

cd  $HOME_DIR
grep ">" $SAMPLE |cut -f 2 -d "_"  |cut -f 1 -d ">" |cut -f 1 -d '/'>Genes.txt
cat Genes.txt |while read line; do rm -rf $line; mkdir $line; done;
Rscript AnaR1.r $HOME_DIR $REF $SAMPLE;















################################################
### Step 3


cd $HOME_DIR;
cat Genes.txt |while read line; do cd $line; mafft --auto $line"_o.fa" >$line"_aln.fa"; cd ..; done;










################################################
### Step 4

#Evaluate the alignment process, check the aligned fast files



################################################
### Step 5

Rscript AnaR2.r $HOME_DIR $SAMPLE_NAME;












################################################
### Step 6

#Refine Alignment


################################################
### Step 7

cd $HOME_DIR;
Rscript AnaR3.r $HOME_DIR $SAMPLE_NAME;











################################################
### Step 8

cd $HOME_DIR;
Rscript AnaR4.r $HOME_DIR $SAMPLE_NAME;










################################################
### Step 9

cd $HOME_DIR;
Rscript AnaR5.r $HOME_DIR $SAMPLE_NAME;










################################################
### Step 10

cd $HOME_DIR;
mkdir CONCAT;
mv ./CONCATINATED_o.fa ./CONCAT/;
mv ./CONCATINATED_aln_e.phy ./CONCAT/;
cd ./CONCAT;







################################################
### Step 11

cd $HOME_DIR;
mkdir TREE_FILES;
cat Genes.txt |while read line; do 
cd $line;
RAND=$(( $RANDOM %99999)); #Random seed
phyml -i $line"_aln_e.phy" -d aa -b 100 -m JTT -c 4 -a e -s BEST -v e -o tlr -f m --rand_start --n_rand_starts 3 --r_seed $RAND --print_site_lnl --print_trace --no_memory_check
cp ./$line"_aln_e.phy_phyml_tree.txt" ../TREE_FILES/;
cd ..;
done;

#This one is for the final concatenated alignment!
RAND=$(( $RANDOM %999));

cd $HOME_DIR;
cd ./CONCAT;
phyml -i >CONCATINATED_aln_e.phy -d aa -b 100 -m JTT -c 4 -a e -s BEST -v e -o tlr -f m --rand_start --n_rand_starts 3 --r_seed $RAND --print_site_lnl --print_trace --no_memory_check;
cp ./CONCATINATED_aln_e.phy_phyml_tree.txt ../TREE_FILES/;

