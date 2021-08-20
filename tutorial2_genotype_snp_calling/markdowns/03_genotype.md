
#### 3. Genotype calling

![stage2](../files/stage2.png)

Here we will explore several ways to call genotypes from sequencing data.
We will also calculate genotypes probabilities to each site for each individual.

In ANGSD, the option to call genotypes is `-doGeno`:
```
angsd -doGeno
...
-doGeno 0
        1: write major and minor
        2: write the called genotype encoded as -1,0,1,2, -1=not called
        4: write the called genotype directly: eg AA,AC etc
        8: write the posterior probability of all possible genotypes
        16: write the posterior probability of called genotype
        32: write the posterior probabilities of the 3 gentypes as binary
        -> A combination of the above can be choosen by summing the values, EG write 0,1,2 types with majorminor as -doGeno 3
        -postCutoff=0.333333 (Only genotype to missing if below this threshold)
        -geno_minDepth=-1       (-1 indicates no cutof)
        -geno_maxDepth=-1       (-1 indicates no cutof)
        -geno_minMM=-1.000000   (minimum fraction af major-minor bases)
        -minInd=0       (only keep sites if you call genotypes from this number of individuals)

        NB When writing the posterior the -postCutoff is not used
        NB geno_minDepth requires -doCounts
        NB geno_maxDepth requires -doCounts
```

Therefore, if we set `-doGeno 2`, genotypes are coded as 0,1,2, as the number of alternate/minor alleles.
If we want to print the major and minor alleles as well then we set `-doGeno 3`.

To calculate the posterior probability of genotypes we need to define a model.
```
angsd -doPost
...
-doPost 0       (Calculate posterior prob 3xgprob)
        1: Using frequency as prior
        2: Using uniform prior
        3: Using SFS as prior (still in development)
        4: Using reference panel as prior (still in development), requires a site file with chr pos major minor af ac an
...
```
`-doPost 2` uses a uniform prior.

Furthermore, this calculation requires the specification of how to assign the major and minor alleles (if biallelic).
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

A typical command for genotype calling is (assuming we analyse our PANY samples):
```
angsd -b $DIR/PANY_bams.txt -ref $REF -out Results/PANY \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
    -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 \
	-GL 2 -doGlf 1

angsd -glf Results/PANY.glf.gz -fai $REF.fai -nInd 15 -out Results/PANY \
	-doMajorMinor 1 -doGeno 3 -doPost 2 -doMaf 1
```
Let's ignore the `-doMaf` option now. We will discuss it later.
This command should take 1 minute to run.

Have a look at the output file:
```
less -S Results/PANY.geno.gz
```
The columns are: chromosome, position, major allele, minor allele, genotypes is 0,1,2 format.

**QUESTION**
How many sites have at least one missing genotype?
```
zcat Results/PANY.geno.gz | grep -1 - | wc -l
```
How many sites do we have?
```
zcat Results/PANY.geno.gz | | wc -l
```
Why is that?

You can control how to set missing genotype when their confidence is low with `-postCutoff`.
For instance, we can set as missing genotypes when their (highest) genotype posterior probability is below 0.95:
```
angsd -glf Results/PANY.glf.gz -fai $REF.fai -nInd 15 -out Results/PANY \
        -doMajorMinor 1 -doGeno 3 -doPost 2 -doMaf 1 -postCutoff 0.95
```

How many sites do we have in total?
How many sites have at least one missing genotype now?
```
zcat Results/PANY.geno.gz | wc -l
zcat Results/PANY.geno.gz | grep -1 - | wc -l
```

Why are there some many sites with missing genotypes?

The mean depth per sample is around 1-2X, therefore genotypes cannot be assigned with very high confidence.
Sites where all genotypes are missing are skipped in the output file.

Setting this threshold depends on the mean sequencing depth of your data, as well as your application.
For some analyses you need to work only with high quality genotypes (e.g. measure of proportion of shared SNPs for gene flow estimate), while for others you can be more relaxed (e.g. estimate of overall nucleotide diversity).
We will show later how to accurately estimate summary statistics with low-depth data.

**BONUS QUESTION**
Try to set the threshold on -postCutoff to 0.50. How many sites do you retrieve? Why? Which value would you choose?


![stage2A](../files/stage2A.png)

--------------------------------

**EXERCISE**

If we assume HWE, then we can use this information as prior probability to calculate genotype posterior probabilities.
The command line would be:
```
angsd -glf Results/PANY.glf.gz -fai $REF.fai -nInd 15 -out Results/PANY \
        -doMajorMinor 1 -doGeno 3 -doPost 1 -doMaf 1
```
using the option `-doPost 1`.

In ANGSD we can restrict our analyses on a subset of positions of interest using the `-sites` option.
The file with these positions need to be formatted as (chromosome positions).
```
echo Mme_chr24:2558528-4558528 48 > Data/snp.txt
echo Mme_chr24:2558528-4558528 61 >> Data/snp.txt
```
We need to index this file in order for ANGSD to process it.
```
angsd sites index Data/snp.txt
```

We are interested in calculating the derived allele frequencies, so are using the ancestral sequence to polarise the alleles.
We also want to compute the allele frequencies for each population separately.
We need to use a different file for each population, with a different list of BAM files, as provided:
```
ls -d -1 $DIR/*_bams.txt
```
```
/home/ubuntu/Share/data/ALL_bams.txt
/home/ubuntu/Share/data/JIGA_bams.txt
/home/ubuntu/Share/data/MAQU_bams.txt
/home/ubuntu/Share/data/MBNS_bams.txt
/home/ubuntu/Share/data/PANY_bams.txt
```
We retain only these populations: Jekyll Island (JIGA), Patchogue (PANY), Minas Basin (MBNS), Magdalen Island (MBNS).

We are calculating the **derived** allele frequencies based on assigned genotypes for each population at these two positions.
We have to specify a putative ancestral sequence, with the variable `$ANC`.

Write the code that performs the following genotype calling for our variants of interest in all populations.
Also, you can directly call genotypes without generating the genotype likelihood files, by starting from bam files directly.
As an indication, you can follow these guidelines:
- use the SAMtools genotype likelihood model
- calculate genotype posterior probabilities using a HWE-based prior
- filter out bases with a quality score less than 20
- filter our reads with a mapping quality score less than 20
- use ony sites where you have at least five samples with data (-mindInd)
- do not set any filtering based on min and max depth
- use -doMajorMinor 5 (try to understand why) and -doMaf 1 options
- set genotypes as missing if the highest genotype probability is less than 0.50
- use option `-sites Data/snp.txt` to restrict the analysis only on selected sites
but feel free to choose some parameters yourself.

```
...
```

Once done, open the output files and calculate the derived allele frequency by counting genotypes.
What is the derived allele frequency for each population for each site?

Can you comment these results?
Do you see any allele frequency differentiation in the derived state by counting genotypes?

If you don't obtain an output, it means that for these loci there is no data. Try to run the calculation on the whole data set and pick one site with data for counting genotypes, if you wish so.

----------------------------

You have now learnt how to call genotypes with ANGSD and appreciated the effect of data uncertainty.
You are now going to learn how to perform SNP calling and estimation of allele frequencies.

[click here](https://github.com/nt246/lcwgs-guide-tutorial/blob/main/tutorial2_genotype_snp_calling/markdowns/04_snp.md) to move to the next session.

------------------------------------



