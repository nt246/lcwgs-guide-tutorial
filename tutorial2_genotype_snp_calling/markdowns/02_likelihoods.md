
#### 2. Genotype likelihoods

![stage1](../files/stage1.png)

We now wish to calculate the ***genotype likelihoods*** for each site at each individual.

To do so you need to specify which genotype likelihood model to use.
```
angsd -GL
...
-GL=0: 
	1: SAMtools
	2: GATK
	3: SOAPsnp
	4: SYK
	5: phys
	6: Super simple sample an allele type GL. (1.0,0.5,0.0)
	7: outgroup gls
	-trim		0		(zero means no trimming)
	-tmpdir		angsd_tmpdir/	(used by SOAPsnp)
	-errors		(null)		(used by SYK)
	-minInd		0		(0 indicates no filtering)

Filedumping:
	-doGlf	0
	1: binary glf (10 log likes)	.glf.gz
	2: beagle likelihood file	.beagle.gz
	3: binary 3 times likelihood	.glf.gz
	4: text version (10 log likes)	.glf.gz
```
A description of these different implementation can be found [here](http://www.popgen.dk/angsd/index.php/Genotype_likelihoods).
The GATK model refers to the first GATK paper, SAMtools is somehow more sophisticated (non-independence of errors), SOAPsnp requires a reference sequence for recalibration of quality scores, SYK is error-type specific.
For most applications and data, GATK and SAMtools models should give similar results.

Let's assume to work with PANY samples only.
A possible command line to calculate genotype likelihoods might be:
```
angsd -b $DIR/PANY_bams.txt -ref $REF -out Results/PANY \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 \
        -GL 2 -doGlf 4
```

where we specify:
* -GL 2: genotype likelihood model as in GATK
* -doGlf 4: output in text format

Ignore the various warning messages. If it is too slow, the add `-nThreads 10` at the end of the command line. This command should take around 2 minutes to run.

![stage1A](../files/stage1A.png)

**QUESTION**
What are the output files?
What's the information inside them?

```
ls Results/PANY.*
```

```
less -S Results/PANY.arg
less -S Results/PANY.glf.gz
```

**BONUS QUESTION**
Try to output files in binary format. Which option should you use? Can you open these files?
Look at the file sizes of text vs binary format. Which one is smaller? Which one would you use?

**BONUS QUESTION**
Try to change some filtering options and record the number of entries in the final output file.


You have learnt how to calculate and read genotype likelihood files.
Now you are going to learn how to perform genotype calling with ANGSD.

[click here](https://github.com/nt246/lcwgs-guide-tutorial/blob/main/tutorial2_genotype_snp_calling/markdowns/03_genotype.md) to move to the next session.

--------------------------------



