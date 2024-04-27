# Step-by-step analysis of paleoproteomic data

## Contents

- [Overview](#overview)
- [Repository Contents](#repository-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Instructions for use](#instructions-for-use)
- [Citation](#citation)

# Overview

This repository contains the scripts and example data for the analysis of enamel proteome sequences from ancient specimens. These scripts can be used to:

1. obtain protein sequences for the enamel proteins starting from sequencing reads mapped to the human reference genome (BAM files)
2. add an ancient protein sequence to a precomputed protein alignment
3. estimate missingness and informativeness for the ancient protein sequence in such alignment
4. mask ambiguous I/L in the protein alignment at relevant positions
5. obtain a simple phylogenetic tree

# Repository Contents

- [Data](./Data): example data.
- [Genes_human](./Genes_human): fasta sequences of human enamel proteins. 
- [Scripts](./Scripts): `R` scripts. 


# System Requirements

- UNIX OS
- R (tested on version 4.0.3 (2020-10-10))
- R packages: "ShortRead", "stringr" , "data.table", "phyclust"
- mafft v.7.205
- PhyML v.3.1 
- ANGSDv0.921 (or newest version)
- BLAST+ v2.7.1 


# Installation Guide

```
git clone https://github.com/johnpatramanis/Nature_Prot_Enamel.git
```

# Demo

**Example on how to use these scripts to reconstruct the ENAM protein sequence for the Altai Neanderthal starting from a BAM file**

1. Define the working directories and software:

Directories:
```
BAM_NEAN='AltaiNea.chr4.bam'
ENAM_HUMAN='Data/ENAM_HUMAN.fa'
STARTS='Data/starts.txt'
GENE_LOCS='Data/Gene_locs.txt'
```

Software:
```
ANGSD='/path/to/angsd'
MAKEBLASTDB='/path/to/makeblastdb'
BLASTALL='/path/to/blastall'
```
2.Create a majority count consensus sequence from the BAM file for the region spanning the enamelin gene 

```
$ANGSD -minQ 20 -minMapQ 30 -doFasta 2 -doCounts 1 -basesPerLine 60 -i $BAM_NEAN -r 4:71494461-71552533 -out Neanderthal_ENAM
```

3. *In silico* splicing of the enamelin gene to remove the introns

```
Rscript Scripts/DataR1.r $GENE_LOCS Neanderthal_ENAM.fa.gz ENAM  Neanderthal

Rscript Scripts/DataR2.r $STARTS Data/ENAM_ei.txt ENAM Neanderthal  
```

4. Format the spliced enamelin gene as a BLAST database

```
$MAKEBLASTDB -dbtype nucl -in Neanderthal_ENAM_spliced.fa
```

5. Perform a BLAST search of the the human ENAM protein against the spliced enamelin gene

```
$BLASTALL  -p tblastn -i $ENAM_HUMAN -d Neanderthal_ENAM_spliced.fa -o Neanderthal_ENAM_spliced.blast -F F -E 32767 -G 32767 -n T -m 0 -M PAM70
```

6. Extract the translated ENAM protein for the Altai Neanderthal from the BLAST result

```
Rscript Scripts/DataR3.r Neanderthal_ENAM_spliced.blast ENAM Neanderthal   
```

# Instructions for use

The step-by-step process for the usage of the scripts is fully described in Taurozzi A. et al 2024.

# Citation

For usage of these scripts cite the following manuscript:

Deep-time phylogenetic inference by palaeoproteomic analysis of dental enamel. Taurozzi A. et al 2024.








