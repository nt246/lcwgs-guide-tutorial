library(tidyverse)

cov <- as.matrix(read.table("/workdir/arne/physalia_lcwgs_data/data_practicals/Results/MME_ANGSD_PCA.covMat", header = F))

pop <- c("JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA"
         ,"PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY"
         ,"MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS"
         ,"MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU")


mme.pca <- eigen(cov)
eigenvectors <- (mme.pca$vectors)
pca.vectors <- as_tibble(cbind(pop, eigenvectors))
x = type_convert(pca.vectors)

pca = ggplot(data = x, aes(x=V2, y=V3, fill = pop, colour = pop)) +
  geom_point(size = 5, shape = 21) +
  xlab("Principal component 1") +
  ylab("Principal component 2")

ggsave(filename = "/workdir/arne/physalia_lcwgs_data/data_practicals/Results/pca_plot.pdf", plot = pca)