
#### 4. SNP calling and allele frequencies

![stage3](../files/stage3.png)

We now want to estimate allele frequencies at each site without relying on genotype calls.
In other words, at each site we want to to estimate (or count) how many copies of different alleles (two in case of biallelic variants) we observe in our sample (across all sequenced individuals).
However with low depth data direct counting of individually assigned genotypes can lead to biased allele frequencies.

ANGSD has an option to estimate **allele frequencies** taking into account data uncertainty from genotype likelihoods:
```
angsd -doMaf
...
-doMaf  0 (Calculate persite frequencies '.mafs.gz')
        1: Frequency (fixed major and minor)
        2: Frequency (fixed major unknown minor)
        4: Frequency from genotype probabilities
        8: AlleleCounts based method (known major minor)
        NB. Filedumping is supressed if value is negative
-doPost 0       (Calculate posterior prob 3xgprob)
        1: Using frequency as prior
        2: Using uniform prior
        3: Using SFS as prior (still in development)
        4: Using reference panel as prior (still in development), requires a site file with chr pos major minor af ac an
Filters:
        -minMaf         -1.000000       (Remove sites with MAF below)
        -SNP_pval       1.000000        (Remove sites with a pvalue larger)
        -rmTriallelic   0.000000        (Remove sites with a pvalue lower)
Extras:
        -ref    (null)  (Filename for fasta reference)
        -anc    (null)  (Filename for fasta ancestral)
        -eps    0.001000 [Only used for -doMaf &8]
        -beagleProb     0 (Dump beagle style postprobs)
        -indFname       (null) (file containing individual inbreedcoeficients)
        -underFlowProtect       0 (file containing individual inbreedcoeficients)
NB These frequency estimators requires major/minor -doMajorMinor
```

Therefore, the estimation of allele frequencies requires the specification of how to assign the major and minor alleles (if biallelic).
```
angsd -doMajorMinor
...
        -doMajorMinor   0
        1: Infer major and minor from GL
        2: Infer major and minor from allele counts
        3: use major and minor from a file (requires -sites file.txt)
        4: Use reference allele as major (requires -ref)
        5: Use ancestral allele as major (requires -anc)
        -rmTrans: remove transitions 0
        -skipTriallelic 0
```

A possible command line to estimate allele frequencies might be (this may take 1 min to run):
```
angsd -b $DIR/PANY_bams.txt -ref $REF -out Results/PANY \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 \
        -GL 1 -doGlf 1 -doMajorMinor 1 -doMaf 1
```
where we specify:
* -doMajorMinor 1: both alleles are inferred from genotype likelihoods
* -doMaf 1: major and minor are fixed

**QUESTION**
What are the output files?

* "Results/PANY.arg"
* "Results/PANY.mafs.gz"

`.args` file is a summary of all options used, while `.mafs.gz` file shows the allele frequencies computed at each site.

Have a look at this file which contains estimates of allele frequencies.
```
zcat Results/PANY.mafs.gz | head
```
and you may see something like
```
chromo	position	major	minor	ref	knownEM	nInd
Mme_chr24:2558528-4558528	27	A	C	A	0.000004	5
Mme_chr24:2558528-4558528	28	T	A	T	0.000004	5
Mme_chr24:2558528-4558528	29	T	A	T	0.000004	5
Mme_chr24:2558528-4558528	30	T	A	T	0.000004	5
Mme_chr24:2558528-4558528	31	A	C	A	0.000004	5
Mme_chr24:2558528-4558528	32	T	A	T	0.000004	5
Mme_chr24:2558528-4558528	35	C	A	C	0.000004	5
Mme_chr24:2558528-4558528	38	A	C	A	0.000004	5
Mme_chr24:2558528-4558528	39	T	A	T	0.000004	6
```

where `knownEM` specifies the algorithm used to estimate the allele frequency which is given under that column.
Please note that this refers to the allele frequency of the allele labelled as `minor`.
The columns are: chromosome, position, major allele, minor allele, reference allele, allele frequency, p-value for SNP calling (if -SNP-pval was called), number of individuals with data.
The last column gives the number of samples with data (you can see that this never below 5 given our filtering).

You can notice that many sites have low allele frequency, probably reflecting the fact that that site is monomorphic.
We may be interested in looking at allele frequencies only for sites that are actually variable in our sample.
Therefore we want to perform a **SNP calling**.

There are two main ways to call SNPs using ANGSD with these options:
```
        -minMaf         0.000000        (Remove sites with MAF below)
        -SNP_pval       1.000000        (Remove sites with a pvalue larger)
```
Therefore we can consider assigning as SNPs sites whose estimated allele frequency is above a certain threhsold (e.g. the frequency of a singleton) or whose probability of being variable is above a specified value.

**QUICK EXERCISE**

As an illustration, call SNPs by computing:
 - genotype likelihoods using GATK method;
 - major and minor alleles inferred from genotype likelihoods;
 - frequency from known major allele but unknown minor;
 - SNPs as those having MAF=>0.05.

Try to write down this command by yourself and comment the results.

```
...

```

As a general guidance, `-GL 1`, `-doMaf 1/2` and `-doMajorMinor 1` should be the preferred choice when data uncertainty is high.
If interested in analysing very low frequency SNPs, then `-doMaf 2` should be selected.
When accurate information on reference sequence or outgroup are available, one can use `-doMajorMinor` to 4 or 5.
Also, detecting variable sites based on their probability of being SNPs is generally a better choice than defining a threshold on the allele frequency.
However, various cutoffs and a dedicated filtering should be perform to assess robustness of your called SNPs.

![stage3A](../files/stage3A.png)

**QUICK EXERCISE**

Try varying the cutoff for SNP calling and record how many sites are predicted to be variable for each scenario.
Identify which sites are not predicted to be variable anymore with a more stringent cutoff (e.g. between a pair of scenario), and plot their allele frequencies.
Use the previously calculated genotype likelihoods as input file (add proper values for ```-glf ? -fai ? -nInd ?```).
```
# iterate over some cutoffs (you can change these)
for PV in 0.05 1e-2 1e-4 1e-6
do
        if [ $PV == 0.05 ]; then echo SNP_pval NR_SNPs; fi
        angsd -glf ? -nInd 15 -fai ? -out Results/PANY.$PV \
                -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
                -SNP_pval $PV &> /dev/null
        echo $PV `zcat Results/PANY.$PV.mafs.gz | tail -n+2 | wc -l`
done
```

A possible output is (your numbers may be different):
```
SNP_pval NR_SNPs
0.05 23052
1e-2 16630
1e-4 9198
1e-6 6348
```

Which sites differ from 0.05 and 0.01? What is their frequency?
This script will also print out the first 20 discordant sites (pK.EM is the p-value for the SNP calling test).
```
Rscript -e 'mafs1 <- read.table(gzfile("Results/PANY.1e-2.mafs.gz"), he=T, row.names=NULL, strings=F); mafs5 <- read.table(gzfile("Results/PANY.0.05.mafs.gz"), header=T, row.names=NULL, stringsAsFact=F); mafs5[!(mafs5[,2] %in% mafs1[,2]),][1:20,]; pdf(file="Results/diffSnpCall.pdf"); par(mfrow=c(1,2)); hist(as.numeric(mafs5[!(mafs5[,2] %in% mafs1[,2]),][,6]), main="Discordant SNPs", xlab="MAF", xlim=c(0,0.5)); hist(as.numeric(mafs5[(mafs5[,2] %in% mafs1[,2]),][,6]), main="Concordant SNPs", xlab="MAF", xlim=c(0,0.5)); dev.off();'
```

You can scp the pdf file to your local machine (from your local machine: scp -i XYZ.pem userX@IP:~/day2/Results/diffSnpCall.pdf .) and visualise it with
```
evince Results/diffSnpCall.pdf
```

What can you conclude from these results?
Which frequencies are more difficult to estimate and therefore affect SNP calling?

----------------------------------

**EXERCISE**

Estimate derived allele frequencies for all populations of interest using a likelihood approach, without relying on genotype calls.
What is the difference compared to what previously estimated by counting genotypes?

----------------------------------

You are now able to calculate genotype likelihoods and allele frequencies and perform genotype and SNP calling with ANGSD.




