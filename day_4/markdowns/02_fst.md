------------------------------------

#### 2. Population genetic differentiation

Here we are going to calculate **allele frequency differentiation** using the FST and PBS (population branch statistic) metric.
Again, we can achieve this by avoid genotype calling using ANGSD.
From the sample allele frequencies likelihoods (.saf files) we can estimate PBS using the following pipeline.

Note that here we use the previously calculated SFS as prior information.
Also, JIGA is our target population, while MAQU and MBNS are reference populations.
If not already done, you should calculate .saf.idx files for each population, as explained in the section above.

The 2D-SFS will be used as prior information for the joint allele frequency probabilities at each site.
From these probabilities we will calculate the FST and population branch statistic (PBS) using JIGA as target population and MAQU and MBNS as reference populations.
Our goal is to detect selection in JIGA in terms of allele frequency differentiation compared to MAQY and MBNS. 

![stats2_bis](../files/stats2_bis.png)

Specifically, we are computing a slinding windows scan, with windows of 10kbp and a step of 1kbp.
This can be achieved using the following commands.

1) This command will compute per-site FST indexes (please note the order of files). The `-whichFst 1` option is preferred Fst estimator for small sample sizes, as it is the case for our test dataset.
```
realSFS fst index Results/MAQU.saf.idx Results/MBNS.saf.idx Results/JIGA.saf.idx -sfs Results/MAQU.MBNS.sfs -sfs Results/MAQU.JIGA.sfs -sfs Results/MBNS.JIGA.sfs -fstout Results/JIGA.pbs -whichFst 1
```
and you can have a look at their values:
```
realSFS fst print Results/JIGA.pbs.fst.idx | less -S
```
where columns are: chromosome, position, (a), (a+b) values for the three FST comparisons, where FST is defined as a/(a+b).
Note that FST on multiple SNPs is calculated as sum(a)/sum(a+b).

![stats2_tris](../files/stats2_tris.png)

2) The next command will perform a sliding-window analysis, where we define the window size using the `-win` option and the step size using the `-step` option. (If you want non-sliding windows then you set `-step` equal to `-win`. Here we are choosing a window size of 10kb with a step size of 1kb.
```
realSFS fst stats2 Results/JIGA.pbs.fst.idx -win 10000 -step 1000 > Results/JIGA.pbs.fst.txt
```

Have a look at the output file:
```
less -S Results/JIGA.pbs.fst.txt
```
The header is:
```
region  chr     midPos  Nsites  Fst01   Fst02   Fst12   PBS0    PBS1    PBS2
```


<br>

Where are interested in the column `PBS2` which gives the PBS values assuming our JIGA population (coded here as 2) being the target population.
Note that negative PBS and FST values are equivalent to 0.

We are also provided with the individual FST values.
You can see that high values of PBS2 are indeed associated with high values of both Fst02 and Fst12 but not Fst01.

We can plot the results along our region of interest, including our inversion breakpoint annotation 
```
R # open R

pbs = read.table("Results/JIGA.pbs.fst.txt", header = T)
pbs$PBS1[which(pbs$PBS1<0)]=0
pbs$PBS2[which(pbs$PBS2<0)]=0
pbs$PBS0[which(pbs$PBS0<0)]=0

head(pbs)

pdf("Results/pbs_plot.pdf")
plot(pbs$midPos, pbs$PBS2, cex=1)
abline(v=1000000, col="red", lwd=3, lty=2)
dev.off() 

q() # close R
```

It will also print out the maximum PBS value observed as this value will be used in the next part.
This script will also plot the PBS variation in JIGA as target population.
Copy the file to your local machine and open it (```scp -i XXX.pem XXX@IP:~/day4/Results/pbs_plot.pdf .```)

**QUESTION**
In which part of our test sequence does JIGA display signatures of selection?

Now you have learnt how estimate FST and PBS from low-coverage data without even assigning allele frequencies.
Next you will learn how to estimate nucleotide diversity 

[click here](https://github.com/nt246/physalia-lcwgs/blob/main/day_4/markdowns/03_diversity.md) to move to the next session.

-------------------------



