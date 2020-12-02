#!/bin/bash

CHR=$1
START_POS=$2
END_POS=$3

if [[ $START_POS -ge $END_POS ]]; then
    echo "ERROR: start position must be smaller than end position." >&2
    exit -1
fi

# Extract sites from input file
TMP_FILE=`mktemp`
awk -v chr=$CHR -v min=$START_POS -v max=$END_POS 'BEGIN{print "snp1\tsnp2\tdist\tr2p\tD\tDp\tr2"} {split($1,pos1,":"); split($2,pos2,":")} pos1[1]==chr && pos1[2]>=min && pos1[2]<=max && pos2[1]==chr && pos2[2]>=min && pos2[2]<=max' | cut -f 1-7 > $TMP_FILE
N=$((`cat $TMP_FILE | wc -l`-1))
if [[ $N -le 0 ]]; then
    echo "ERROR: no SNPs found in region." >&2
    exit -2
else
    echo "`bc <<< "(1+sqrt(1-4*2*-$N))/2"` SNPs found!" >&2
fi

# Plot LD blocks
cat <<EOF | R --vanilla --slave
library(gtools)
library(reshape2)
library(LDheatmap)

r <- read.table("$TMP_FILE", header=TRUE, stringsAsFactors=FALSE)
id <- unique(mixedsort(c(r[,"snp1"],r[,"snp2"])))
posStart <- head(id,1)
posEnd <- tail(id,1)
r <- rbind(r, c(posStart,posStart,0,NA,NA,NA,NA), c(posEnd,posEnd,0,NA,NA,NA,NA))

for (ld in c("r2")) {
  m <- apply(acast(r, snp1 ~ snp2, value.var=ld, drop=FALSE),2,as.numeric)
  rownames(m) <- colnames(m)
  m <- m[mixedorder(rownames(m)),mixedorder(colnames(m))]
  id <- rownames(m)
  dist <- as.numeric(sub(".*:","",id))

  # Save plot
  pdf(paste("LD_blocks", ld,"pdf", sep="."), width=10, height=10)
  LDheatmap(m, genetic.distances=dist, geneMapLabelX=0.75, geneMapLabelY=0.25, color="blueToRed", LDmeasure=ld)
  x <- dev.off()
}
EOF
