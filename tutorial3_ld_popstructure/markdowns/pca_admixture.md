Tutorial 3b: Population structure and admixture analysis
================

  - [Today’s data](#todays-data)
  - [Initial preparation](#initial-preparation)
  - [Principal components analysis
    (PCA)](#principal-components-analysis-pca)
      - [PCA with LD-pruned SNP
        dataset](#pca-with-ld-pruned-snp-dataset)
          - [Question](#question)
      - [2. Alternative approach: Covariance matrix estimation with
        PCAngsd](#2-alternative-approach-covariance-matrix-estimation-with-pcangsd)
          - [Question](#question-1)
      - [Optional: Compare the results to a PCA based on all
        SNPs](#optional-compare-the-results-to-a-pca-based-on-all-snps)
          - [Question:](#question-2)
          - [Question:](#question-3)
      - [Admixture analysis](#admixture-analysis)
          - [ngsAdmix](#ngsadmix)
              - [Question:](#question-4)
          - [Alternative approach:
            PCAngsd](#alternative-approach-pcangsd)

<br> <br>

In this session you will learn how to use low-coverage whole genome data
to do:

  - Principal Components Analysis (PCA)
  - Admixture analysis

<br> <br>

# Today’s data

As outlined in the previous exercises, we are working with low-coverage
NGS data for 60 Atlantic silversides from the following 4 populations:

<img src="https://github.com/nt246/physalia-lcwgs/blob/main/day_3/img/Silverside_Sample_Map.png?raw=true" height="400">

These populations have been previously studied in [Therkildsen et
al. 2019](https://science.sciencemag.org/content/365/6452/487) and
[Wilder et
al. 2020](https://onlinelibrary.wiley.com/doi/10.1002/evl3.189), and
cover the entire distribution range of the species.

Our NGS data are in bam format and span a 2Mb region on chromosome 24.
The interesting aspect about chromosome 24 is that it harbours a large
polymorphic inversion that differs in its frequency across populations.
The test dataset spans one breakpoint of this inversion (1Mb up and
downstream).

In this tutorial, we want to use these data to examine the population
structure of Atlantic silverside along its distribution range using
principal components analysis (PCA) and admixture analysis. In addition,
we will use principal components analysis and ancestry analyses to
assess which populations most likely contain alternative inversion
karyotypes.

For this practical, we will be using
[ANGSD](http://popgen.dk/wiki/index.php/ANGSD) (Analysis of Next
Generation Sequencing Data),
[ngsAdmix](http://www.popgen.dk/software/index.php/NgsAdmix) and
[PCAngsd](http://www.popgen.dk/software/index.php/PCAngsd).

<br> <br>

# Initial preparation

Again, make sure you have cloned this GitHub repository to your Linux
server, ideally to your home directory, and change your current working
directory to the `tutorial3_ld_popstructure` folder. Set this path as
your `BASEDIR`. For example:

``` bash

BASEDIR=~/lcwgs-guide-tutorial/tutorial3_ld_popstructure/
cd $BASEDIR
ls
```

Then, create a new folder named `results` in the `BASEDIR` for storing
our results.

``` bash

mkdir results
```

And then we set all the rest of the paths:

``` bash

REF=$BASEDIR/reference/Ref.fa
ANGSD=/programs/angsd0.930/angsd/angsd
NGSADMIX=/programs/NGSadmix/NGSadmix ## This is the path to NGSadmix on Cornell BioHPC servers. Make sure that you change this when running on a different server.
SAMTOOLS=samtools ## This is the path to samtools on Cornell BioHPC servers. Make sure that you change this when running on a different server.
PCANGSD=/programs/pcangsd-0.98/pcangsd.py
```

We will also need to index the reference file using samtools

``` bash

$SAMTOOLS faidx $REF
```

<br> <br>

# Principal components analysis (PCA)

To perform PCA with low-coverage NGS data, we have to infer the genetic
covariance matrix, which can be estimated in different ways. Here we
will estimate the covariance matrix using single-read sampling in
[ANGSD](http://www.popgen.dk/angsd/index.php/PCA_MDS) as discussed
during the lecture, but will also show you how to estimate the
covariance matrix using
[PCAngsd](http://www.popgen.dk/software/index.php/PCAngsd)

The first question is what our input dataset should look like?

We only want to focus on variant sites in our population structure
analyses. As we learned yesterday, we could just perform SNP calling in
ANGSD and get a list of all variant sites in our dataset and provide
ANGSD this list of variant sites in future runs using the -sites option.

Here is some example code to illustrate how we would do this. WE ARE NOT
RUNNING THIS CODE TODAY - JUST READ OVER IT, DON’T COPY AND RUN IT.

``` bash

# cd $BASEDIR
# $ANGSD -b $BASEDIR'/sample_lists/ALL_bams.txt' -anc $REF -out $BASEDIR'/results/MME_SNPs' \
# -minMapQ 20 -minQ 20 -doMaf 1 -minMaf 0.05 -SNP_pval 2e-6 \
# -GL 1 -doGlf 2 -doMajorMinor 1 -doPost 1
```

<br>

We could then extract a list of all variant sites from the minor allele
frequency file:

``` bash

# gunzip -c $BASEDIR'/results/MME_SNPs.mafs.gz' | cut -f 1,2,3,4 | tail -n +2 > Global_SNPList_MME_SNPs.txt
```

<br>

And then we have to index our SNP list so that ANGSD can read it:

``` bash

# $ANGSD sites index $BASEDIR'/results/Global_SNPList_MME_SNPs.txt'
```

<br>

However, in general it is good practice to perform LD pruning or at
least thinning to reduce the impact of non-independent SNP clusters,
e.g. large LD clusters in inversions, on the inferred population
structure. In this practical, we will perform the population structure
analysis using an LD-pruned SNP dataset that you prepared earlier today.
Later on we will use the full SNP dataset. The full SNP dataset will
help us to understand how the inversion karyotype is distributed in our
dataset and will highlight the impact of large LD clusters on the
inferred structure.

<br>

## PCA with LD-pruned SNP dataset

**1. Estimating the covariance matrix using a single read sampling
approach:**

Using the list of LD-pruned variant sites and the code shown below, we
can estimate a covariance matrix using single-read sampling for
LD-pruned SNPs using ANGSD.

Let’s specify our SNP list and index it:

``` bash

SNPLIST=$BASEDIR/ngsld/LDpruned_snps.list
$ANGSD sites index $SNPLIST
```

``` bash

$ANGSD -b $BASEDIR'/sample_lists/ALL_bams.txt' -anc $REF -out $BASEDIR'/results/MME_ANGSD_PCA_LDpruned' \
    -GL 1 -doGlf 2 -doMajorMinor 3 -doMAF 1 -doPost 1 -doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1 -sites $SNPLIST
    
```

At the same time, we will also output the genotype likelihoods for these
variant sites in beagle likelihood file format (beagle.gz) by specifying
the `-doGlf 2` option. This will be used as input for estimating the
covariance matrix using PCAngsd (below)

<br>

After angsd finished running, we can see that it produced a range of
different output files:

``` bash

ls results/
```

    MME_ANGSD_PCA_LDpruned.arg  
    MME_ANGSD_PCA_LDpruned.beagle.gz  
    MME_ANGSD_PCA_LDpruned.covMat  
    MME_ANGSD_PCA_LDpruned.ibs.gz  
    MME_ANGSD_PCA_LDpruned.ibsMat  
    MME_ANGSD_PCA_LDpruned.mafs.gz

For the PCA, we are most interested in the `.covMat` file, which
contains the covariance matrix. The other files are an IBS matrix
(`.ibsMat`) that could be used for a multidimensional scaling analysis
(MDS), a file with the minor allele frequencies (`.mafs.gz`) and the
`beagle.gz` genotype likelihood file.

<br>

We can look at the covariance matrix in R, where we will also perform
the PCA using the `eigen` function.

We have to load the covariance matrix into R and then we can optionally
provide population assignments for each individual (rows in same order
as input bam file list `-b $BASEDIR/sample_lists/ALL_bams.txt`). We
recommend that you use R-Studio server for this, but running R on
command line also work (you can start R in the server by typing `R`, and
end the session by typing `quit()`)

``` r
## load the tidyverse package
library(tidyverse)

## define your project directory
basedir <- "~/lcwgs-guide-tutorial/tutorial3_ld_popstructure/" # Make sure to edit this to match your $BASEDIR

#Load the covariance matrix. Don't forget to change the file path if you have downloaded the data to your own computer.
cov <- as.matrix(read.table(paste0(basedir, "results/MME_ANGSD_PCA_LDpruned.covMat"), header = F))

#We will also add a column with population assingments
pop <- c("JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA"
         ,"PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY"
         ,"MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS"
         ,"MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU")

mme.pca <- eigen(cov) #perform the pca using the eigen function. 
```

We can then extract the eigenvectors from the pca object and format them
into a dataframe for plotting, e.g. using `ggplot()`.

``` r
eigenvectors = mme.pca$vectors #extract eigenvectors 
pca.vectors = as_tibble(cbind(pop, data.frame(eigenvectors))) #combine with our population assignments

pca = ggplot(data = pca.vectors, aes(x=X1, y=X2, colour = pop)) + geom_point()

ggsave(filename = paste0(basedir, "/results/pca_LDpruned_plot.pdf"), plot = pca) #change file path if data on your own computer
```

Additionally, we can extract the eigenvalues for each eigenvector, and
can then estimate the variance explained for each eigenvector (e.g. here
for PC1 to PC4):

``` r
pca.eigenval.sum = sum(mme.pca$values) #sum of eigenvalues
varPC1 <- (mme.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (mme.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (mme.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (mme.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4
```

<br>

#### Question

How strongly are the four different populations separated? And how much
variation do the first to PCs explain?

<br>

## 2\. Alternative approach: Covariance matrix estimation with PCAngsd

The covariance matrix can also be inferred from genotype likelihoods
using PCAangsd. PCAngsd takes as input genotype likelihoods in beagle
format, which we generated in the step before using the `-doGLF 2`
option.

PCAngsd provides a multitude of different settings, described
[here](http://www.popgen.dk/software/index.php/PCAngsd). We won’t change
any of the settings here and only use the default settings, which are
sufficient in most cases.

We provide the path to the input file using the `-beagle` option, which
also tells PCAngsd that we are working with a beagle file. The output
path and output name is provided using the `-o` option.

``` bash
python3 $PCANGSD -beagle $BASEDIR/results/MME_ANGSD_PCA_LDpruned.beagle.gz -o $BASEDIR/results/PCAngsd_LDpruned_covmat
```

<br>

**Optional** We can perform the principal components analysis and plot
PC1 vs PC2 the same way we did before.

``` r
#Load the RcppCNPy package
library(RcppCNPy)

#Define basedir
basedir <- "~/lcwgs-guide-tutorial/tutorial3_ld_popstructure/" # Make sure to edit this to match your $BASEDIR
setwd(basedir)
#Load the covariance matrix
cov <- npyLoad("results/PCAngsd_LDpruned_covmat.cov.npy")

#We will also add a column with population assingments
pop <- c("JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA"
         ,"PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY"
         ,"MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS"
         ,"MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU")
     
mme.pca <- eigen(cov) #perform the pca using the eigen function. 

eigenvectors = mme.pca$vectors #extract eigenvectors 
pca.vectors = as_tibble(cbind(pop, data.frame(eigenvectors))) #combine with our population assignments

pca = ggplot(data = pca.vectors, aes(x=X1, y=X2, colour = pop)) + geom_point()

ggsave(filename = paste0(basedir, "/results/pca_LDpruned_pcangsd_plot.pdf"), plot = pca) #change file path if data on your own computer

pca.eigenval.sum = sum(mme.pca$values) #sum of eigenvalues
varPC1 <- (mme.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (mme.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (mme.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (mme.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4
```

<br>

#### Question

Do the inferred population structures differ between the two approaches?
How different are the inferred variances explained by PC1 and PC2?

**DISCLAIMER:** The inferred patterns and differences between the
results will be strongly impacted by the small size of our test region.

-----

<br>

## Optional: Compare the results to a PCA based on all SNPs

While it is best practice to perform a PCA based on an LD-pruned SNP
dataset, PCA are often performed based on full SNP datasets. If you are
interested and have time, you could compare the PCA for LD-pruned SNPs
to one that was performed for all SNPs in the dataset without
LD-pruning.

For this, we will only estimate the covariance matrix using single read
sampling in ANGSD. SNPs can be inferred based on genotype likelihoods
(see day 2) at the same time as inferring the covariance matrix
(`-doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1`) by providing the
`-SNP_pval` option. SNPs should be restricted to more common variants
with minor allele frequencies of at least 5% using the `-minMAF 0.05`
option.

``` bash
$ANGSD -b $BASEDIR/sample_lists/ALL_bams.txt -anc $REF -out $BASEDIR'/results/MME_ANGSD_PCA' \
    -minMapQ 20 -minQ 20 -doMaf 1 -minMaf 0.05 -SNP_pval 2e-6 \
    -GL 1 -doGlf 2 -doMajorMinor 1 -doPost 1 \
    -doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1 -P 4
```

<br>

Again, we perform the principal components analysis using the `eigen`
function in R:

``` r
library(tidyverse)
#Define basedir
basedir <- "~/lcwgs-guide-tutorial/tutorial3_ld_popstructure/" # Make sure to edit this to match your $BASEDIR

#Load the covariance matrix
cov <- as.matrix(read.table(paste0(basedir, "results/MME_ANGSD_PCA.covMat"), header = F))

#We will also add a column with population assingments
pop <- c("JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA"
         ,"PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY"
         ,"MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS"
         ,"MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU")
     
mme.pca <- eigen(cov) #perform the pca using the eigen function. 

eigenvectors = mme.pca$vectors #extract eigenvectors 
pca.vectors = as_tibble(cbind(pop, data.frame(eigenvectors))) #combine with our population assignments

pca = ggplot(data = pca.vectors, aes(x=X1, y=X2, colour = pop)) + geom_point()

ggsave(filename = paste0(basedir, "results/pca_allSNPs_plot.pdf"), plot = pca) #change file path if data on your own computer
```

<br>

#### Question:

How do population relationships differ between the on LD-pruned PCA and
the full-dataset PCA? Why do the population relationships differ?

<br>

<details>

<summary>Click here for a hint</summary>

<br>

As you can see the separation between the populations, particularly
between JIGA and the other populations, is a lot stronger, without much
variation within populations. The reason for that is that JIGA has a
different karyotype for the inversion 24 compared to MAQU, MBNS and
PANY. Variation along PC2 represents genetic variation within each
inversion karyotype.

</details>

<br> <br>

#### Question:

What’s with the lonely PANY individual that is intermediate between the
two main clusters?

<br>

-----

<br>

## Admixture analysis

In some cases we also want to infer genome-wide admixture proportions
for each individuals. Similar to PCA, there are different ways of
inferring admixture proportions from genotype likelihoods. Here, we will
use [ngsAdmix](http://www.popgen.dk/software/index.php/NgsAdmix) and
will also present a way of inferring admixture proportions with PCAngsd.
Similar to the PCA, we want to use an LD-pruned SNP dataset for this
analysis to reduce the impact of non-independent SNPs on the ancestry
inference.

<br>

### ngsAdmix

ngsAdmix uses a genotype likelihood file in beagle format (same as for
PCAngsd) as input, which is specified using the `-likes` option. In
addition, there are a range of parameters that can be adjusted. Here we
only set the number of ancestry clusters using the `-K` option to K=2.
In reality, it is advisable to compare different numbers of ancestry
clusters by iterating over different values of K.

``` bash

$NGSADMIX -likes $BASEDIR'/results/MME_ANGSD_PCA_LDpruned.beagle.gz' -K 2 -o $BASEDIR'/results/MME_LDpruned_ngsAdmix_K2_out'
```

<br>

In case the analysis is not converging, one can also increase the
maximum number of EM iterations using the `-maxiter` option.

ngsAdmix produces three different outputs files:

  - A run log containing the log likelihood of the admixture estimates:
    `.log file`
  - A file containing the allele frequencies of all ancestral
    populations (in this case two ancestral clusters): `.fopt file`
  - A file containing the admixture proportions for each individual:
    `.qopt file`

We are mostly interested in the admixture proportions and the log
likelihood of the estimates. The log likelihoods can be compared between
runs with different values of K to select the most likely number of
ancestral clusters (However, this should always be interpreted in a
biologically meaningful context)

<br>

#### Question:

What is the likelihood for `-K 2`?

<details>

<summary>Click here to expand</summary>

``` bash

cat $BASEDIR/results/MME_LDpruned_ngsAdmix_K2_out.log
```

    Input: lname=results/MME_ANGSD_PCA_LDpruned.beagle.gz nPop=2, fname=(null) qname=(null) outfiles=results/MME_LDpruned_ngsAdmix_K2_out
    Setup: seed=1603211682 nThreads=1 method=1
    Convergence: maxIter=2000 tol=0.000010 tolLike50=0.100000 dymBound=0
    Filters: misTol=0.050000 minMaf=0.050000 minLrt=0.000000 minInd=0
    Input file has dim: nsites=247 nind=60
    Input file has dim (AFTER filtering): nsites=245 nind=60
        [ALL done] cpu-time used =  0.11 sec
        [ALL done] walltime used =  0.00 sec
    best like=-11324.127798 after 93 iterations

</details>

<br>

Furthermore, we can plot admixture proportions in R:

``` r
basedir <- "~/lcwgs-guide-tutorial/tutorial3_ld_popstructure/" # Make sure to edit this to match your $BASEDIR
library(tidyverse) #load the tidyverse package for formatting and plotting
```

    ## Warning: replacing previous import 'lifecycle::last_warnings' by
    ## 'rlang::last_warnings' when loading 'pillar'

    ## Warning: replacing previous import 'lifecycle::last_warnings' by
    ## 'rlang::last_warnings' when loading 'tibble'

    ## Warning: replacing previous import 'lifecycle::last_warnings' by
    ## 'rlang::last_warnings' when loading 'hms'

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✓ ggplot2 3.3.3     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.2     ✓ dplyr   1.0.7
    ## ✓ tidyr   1.1.3     ✓ stringr 1.4.0
    ## ✓ readr   1.4.0     ✓ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
#Load the covariance matrix
admix = read_table(paste0(basedir, "/results/MME_LDpruned_ngsAdmix_K2_out.qopt"), col_names = F)
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   X1 = col_double(),
    ##   X2 = col_double()
    ## )

``` r
#We will also add a column with population assingments
pop <- c("JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA"
         ,"PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY"
         ,"MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS"
         ,"MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU")

admix.id = as.data.frame(cbind(pop, admix))
names(admix.id) = c("pop","q1","q2")

pdf(paste0(basedir, "/results/NGSadmix_LDpruned_K2_plot.pdf"))
plot = barplot(t(as.matrix(subset(admix.id, select=q1:q2))), col=c("firebrick","royalblue"), border=NA)
dev.off() 
```

    ## png 
    ##   2

**Optional** If you want, you can try changing the value for `-K` and
compare the inferred log likelihoods and admixture proportions. Which
value of K has the stronger support?

<br>

``` bash

$NGSADMIX -likes $BASEDIR/results/MME_ANGSD_PCA_LDpruned.beagle.gz -K 3 -o $BASEDIR/results/MME_LDpruned_ngsAdmix_K3_out
```

<br>

<details>

<summary>Click here to expand</summary>

``` bash

cat $BASEDIR/results/MME_LDpruned_ngsAdmix_K3_out.log
```

    Input: lname=results/MME_ANGSD_PCA_LDpruned.beagle.gz nPop=3, fname=(null) qname=(null) outfiles=results/MME_LDpruned_ngsAdmix_K3_out
    Setup: seed=1603212799 nThreads=1 method=1
    Convergence: maxIter=2000 tol=0.000010 tolLike50=0.100000 dymBound=0
    Filters: misTol=0.050000 minMaf=0.050000 minLrt=0.000000 minInd=0
    Input file has dim: nsites=247 nind=60
    Input file has dim (AFTER filtering): nsites=245 nind=60
        [ALL done] cpu-time used =  1.02 sec
        [ALL done] walltime used =  1.00 sec
    best like=-11077.337705 after 750 iterations

</details>

<br>

### Alternative approach: PCAngsd

In addition to ngsAdmix, PCAngsd can also be used to infer admixture
proportions. The command is similar to the one used for the PCA with the
addition of to a `-admix` option. The input file for ngsAdmix and
PCAngsd is the same, making comparisons relatively easy.

Other than ngsAdmix, PCAngsd can automatically infer the most likely
number of ancestry cluster. However, one can also set the number of
clusters using the `-admix_K` option.

<br>

``` bash

python3 $PCANGSD -beagle $BASEDIR/results/MME_ANGSD_PCA_LDpruned.beagle.gz -admix -admix_K 2 -o $BASEDIR/results/MME_PCAngsd_K2_out
```

<br>

**Optional** If you want to, you can plot these results as well and
compare them to the estimates with NGSadmix.
