Tutorial 3a: Linkage disequilibrium estimation
================

  - [Define paths to the project directory and
    programs](#define-paths-to-the-project-directory-and-programs)
      - [Set the project directory as a variable named
        `BASEDIR`](#set-the-project-directory-as-a-variable-named-basedir)
      - [Specify the path to required programs as
        variables](#specify-the-path-to-required-programs-as-variables)
  - [Estimate LD](#estimate-ld)
      - [Prepare the input files](#prepare-the-input-files)
      - [Run ngsLD](#run-ngsld)
      - [Visualize LD blocks](#visualize-ld-blocks)
          - [Install required R packages on your local
            computer](#install-required-r-packages-on-your-local-computer)
          - [Run `LD_blocks.sh`](#run-ld_blockssh)
      - [LD pruning](#ld-pruning)
          - [Run LD pruning](#run-ld-pruning)
          - [Generate an LD-pruned SNP
            list](#generate-an-ld-pruned-snp-list)

<br> <br>

The estimation of linkage disequilibrium (LD) has important
applications, e.g. for inference of population size, demographic
history, selection, and for the discovery of structural variants. In
addition, since many downstream analyses make assumptions about the
independence of genomic loci, LD estimates are essential for trimming
the list of loci to be included in these analyses (LD pruning). In this
session, you will learn to estimate linkage disequilibrium using the
program [ngsLD](https://github.com/fgvieira/ngsLD), which employs a
maximum likelihood approach to account for genotype uncertainty in
low-coverage whole genome sequencing data. You will then visualize the
LD pattern and generate a list of LD-pruned SNPs for downstream
analyses.

We will continue working with low-coverage NGS data for 60 Atlantic
silversides from different four populations studied in [Therkildsen et
al. 2019](https://science.sciencemag.org/content/365/6452/487) and
[Wilder et
al. 2020](https://onlinelibrary.wiley.com/doi/10.1002/evl3.189). The
sequencing data from these individuals have been mapped to a 2Mb section
of chromosome 24. The entire genome is \~650 Mb, so this is just a small
snippet, but it’s an interesting one because it harbours a large
polymorphic inversion that varies substantially in frequency across the
species distribution range. Our test dataset spans one breakpoint of
this inversion (1Mb up and downstream). Therefore, we might expect that
LD patterns dramatically differ across the inversion breakpoint in our
dataset, as inversions can strongly suppress recombination.

<br>

# Define paths to the project directory and programs

We need to make sure the server knows where to find the programs we’ll
be running and our input and output directories. This will always need
to be specified every time we run our scripts in a new login session.

<br>

## Set the project directory as a variable named `BASEDIR`

Again, make sure you have cloned this GitHub repository to your Linux
server, ideally to your home directory, and change your current working
directory to the `tutorial3_ld_popstructure` folder. Set this path as
your `BASEDIR`. For example:

``` bash

BASEDIR=~/lcwgs-guide-tutorial/tutorial3_ld_popstructure
cd $BASEDIR
ls
```

<br>

## Specify the path to required programs as variables

Save the path to ngsLD as a variable named `NGSLD`. For example:

``` bash

NGSLD=/programs/ngsLD/ ## This is the path to ngsLD on Cornell BioHPC servers. Make sure that you change this when running on a different server.
```

Note that this is the path to the folder that contains the `ngsLD`
executable, but not the `ngsLD` executable itself. The reason is that we
will use other scripts with this folder.

<br>

# Estimate LD

## Prepare the input files

ngsLD requires two input files.

1.  `--geno FILE`: input file with genotypes, genotype likelihoods or
    genotype posterior probabilities. With low-coverage data, a genotype
    likelihood file is often used, in which each row is a SNP and each
    individual has three columns correponding to the likelihood of the
    three genotypes (we only need to keep track of the likelihood of
    three genotypes (rather than all ten possible) when the major and
    minor allele has been inferred). Therefore, a `beagle` formatted
    genotype likelihood file generated from ANGSD (`-doGlf 2`) can be
    inputted into ngsLD after the header row and the first three columns
    (i.e. positions, major allele, minor allele) are removed. Here,
    because of time constraint, we will subsample a beagle file
    generated from a similar ANGSD run to what we explored yesterday
    (`MME_ANGSD_PCA.beagle.gz`) as the input to ngsLD by selecting one
    SNP in every 50 SNPs.

2.  `--pos FILE`: input file with site coordinates (one per line), where
    the 1st column stands for the chromosome/contig and the 2nd for the
    position (bp). One convenient way to generate this is by selecting
    the first two columns of the `mafs` file outputted by ANGSD, with
    the header removed. Again, for this exercise, we will downsample the
    `mafs` file that were generated in the same ANGSD run that generated
    our genotype likelihood beagle file (`MME_ANGSD_PCA.mafs.gz`).

<!-- end list -->

``` bash

## Prepare a geno file by subsampling one SNP in every 50 SNPs in the beagle file
zcat $BASEDIR/angsd/MME_ANGSD_PCA.beagle.gz | awk 'NR % 50 == 0' | cut -f 4- | gzip  > $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.beagle.gz

## Prepare a pos file by subsampling one SNP in every 50 SNPs in the mafs filre
zcat $BASEDIR/angsd/MME_ANGSD_PCA.mafs.gz | cut -f 1,2 |  awk 'NR % 50 == 0' | sed 's/:/_/g'| gzip > $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.pos.gz
```

Don’t worry if you don’t understand the `sed` command. The `:` in the
chromosome name interferes with the code, so the `sed` command just
replaces the `:` with `_`.

<br>

## Run ngsLD

Important ngsLD parameters:

  - `--probs`: specification of whether the input is genotype
    probabilities (likelihoods or posteriors)?
  - `--n_ind INT`: sample size (number of individuals).
  - `--n_sites INT`: total number of sites.
  - `--max_kb_dist DOUBLE`: maximum distance between SNPs (in Kb) to
    calculate LD. Set to 0(zero) to disable filter. \[100\]
  - `--max_snp_dist INT`: maximum distance between SNPs (in number of
    SNPs) to calculate LD. Set to 0 (zero) to disable filter. \[0\]
  - `--n_threads INT`: number of threads to use. \[1\]
  - `--out FILE`: output file name. \[stdout\]

<!-- end list -->

``` bash

$NGSLD/ngsLD \
--geno $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.beagle.gz \
--pos $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.pos.gz \
--probs \
--n_ind 60 \
--n_sites 1134 \
--max_kb_dist 0 \
--n_threads 1 \
--out $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.ld 
```

Note: this step may take a few minutes to run. In the meantime, you can
take a look at the [ngsLD GitHub
page](https://github.com/fgvieira/ngsLD).

<br>

## Visualize LD blocks

Now, we’ll visualize the LD estimation results using an R script.

<br>

### Install required R packages on your local computer

Run the following commands in R. You might need to upgrade your base R
to a version \> 4.0 to be able to install these packages.

``` r
install.packages("LDheatmap")
install.packages("reshape2")
install.packages("gtools")
```

<details>

<summary> Click here if you are having trouble installing LDheatmap
</summary>

If you see the following error message when installing `LDheatmap`

    ERROR: dependency ‘snpStats’ is not available for package ‘LDheatmap’

You can try to install `snpStats` using

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("snpStats")
```

Once `snpStats` is installed, you can run
`install.packages("LDheatmap")` again.

</details>

<br>

### Run `LD_blocks.sh`

We will use a script slightly modified from the original `LD_blocks.sh`
script provided by ngsLD to generate a plot of LD blocks in our data. It
takes three argument in the following order:

1.  chromosome / linkage group / scaffold name
2.  starting position
3.  ending position

Run the following script in command line.

``` bash

cd  $BASEDIR/ngsld

cat $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.ld | bash $BASEDIR/scripts/LD_blocks.sh \
Mme_chr24_2558528-4558528 \
200000 \
1400000
```

<br>

Examine the pdf file that the script outputs. If you have not got the
code to work, see the [plot that we have
generated](https://github.com/nt246/lcwgs-guide-tutorial/blob/main/tutorial3_ld_popstructure/ngsld/LD_blocks.r2.pdf).

Do you notice any interesting pattern in this plot of LD blocks? What do
you think is causing this pattern and why?

<br>

## LD pruning

### Run LD pruning

For many downstream analyses, independence among different SNPs is often
assumed, so it is important to generate a list of SNPs that are in low
LD with each other. To do this, we can use the `prune_graph.pl` script
provided in ngsLD. This script takes the LD estimation output (in our
case `MME_ANGSD_PCA.ld`) as its input. Some important parameters
include:

  - `--max_kb_dist INT`: Maximum distance between nodes (ie. SNPs, input
    file 3rd column) to assume they are connected
  - `--min_weight FLOAT`: Minimum weight (in –weight\_field) of an edge
    to assume nodes are connected (the weight refers to the LD estimate
    between two SNPs)
  - `--out FILE`: Path to output file \[STDOUT\]

<br>

``` bash

perl $NGSLD/scripts/prune_graph.pl \
--in_file $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.ld \
--max_kb_dist 2000 \
--min_weight 0.5 \
--out $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled_unlinked.id
```

Check the LD pruning result. How many SNPs survived the LD pruning
process and how many were lost?

<br>

### Generate an LD-pruned SNP list

We will use R to generate an LD-pruned SNP list in a format that can be
used by ANGSD for downstream analyses.

``` r
basedir="~/lcwgs-guide-tutorial/tutorial3_ld_popstructure/"

pruned_position <- as.integer(gsub("Mme_chr24_2558528-4558528:", "", readLines(paste0(basedir, "ngsld/MME_ANGSD_PCA_subsampled_unlinked.id"))))

snp_list <- read.table(paste0(basedir, "angsd/MME_ANGSD_PCA.mafs.gz"), stringsAsFactors = F, header = T)[,1:4]

pruned_snp_list <- snp_list[snp_list$position %in% pruned_position, ]
  
write.table(pruned_snp_list, paste0(basedir, "ngsld/LDpruned_snps.list"), col.names = F, row.names = F, quote = F, sep = "\t")
```
