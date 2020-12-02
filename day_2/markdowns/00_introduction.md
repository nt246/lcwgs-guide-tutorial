
In this session you will learn how to do:
* genotype calling
* allele frequency estimation
* variant (or SNP) calling

We are using the program [ANGSD](http://popgen.dk/wiki/index.php/ANGSD) (Analysis of Next Generation Sequencing Data).
More information about its rationale and implemented methods can be found [here](http://www.ncbi.nlm.nih.gov/pubmed/25420514).

According to its website *ANGSD is a software for analyzing next generation sequencing data. The software can handle a number of different input types from mapped reads to imputed genotype probabilities. Most methods take genotype uncertainty into account instead of basing the analysis on called genotypes. This is especially useful for low and medium depth data.*

Please make sure to follow these preparatory instructions below before running these examples. 
Briefly, you need to set the path to the software and various data that will be used.
Also, you will have to create two folders on your working directory, one for your results and one for your intermediate data.

Make sure you are in your home directory.
```
cd ~
```
and create a folder for this session and enter it
```
mkdir day2
cd day
```
and you should be in `~/day2`.
Also, you will have to create two folders on your working directory, one for your results and one for your intermediate data.
```
mkdir Results
mkdir Data
```

Let's set all paths
```
DIR=/home/ubuntu/Share/data
DATA=$DIR/BAMS
REF=$DIR/Ref.fa
ANC=$DIR/outgrp_ref.fa
```

The **workflow** for this session looks like this

![stages](../files/stages.png)

which seems daunting! 
However, that's not the case and we will go through each step to understand each one of them.

The **workflow** is roughly divided into four steps:

[1.](https://github.com/nt246/physalia-lcwgs/blob/main/day_2/markdowns/01_filtering.md) Data filtering and I/O

[2.](https://github.com/nt246/physalia-lcwgs/blob/main/day_2/markdowns/02_likelihoods.md) Genotype likelihoods

[3.](https://github.com/nt246/physalia-lcwgs/blob/main/day_2/markdowns/03_genotype.md) Genotype calling

[4.](https://github.com/nt246/physalia-lcwgs/blob/main/day_2/markdowns/04_snp.md) SNP calling

You are now going to learn how to build your first pipeline in ANGSD for data processing and filtering.

[click here](https://github.com/nt246/physalia-lcwgs/blob/main/day_2/markdowns/01_filtering.md) to move to the next session.

-----------------------------------------------



