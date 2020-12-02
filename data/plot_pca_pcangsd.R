library(tidyverse) #load the tidyverse package for formatting and plotting

#Load the covariance matrix
cov = as.matrix(RcppCNPy::npyLoad("/Users/arnejacobs/Dropbox/covmatrix.cov.npy"))

#We will also add a column with population assingments
pop <- c("JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA"
         ,"PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY"
         ,"MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS"
         ,"MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU")

mme.pca <- eigen(cov) #perform the pca using the eigen function. 

eigenvectors <- (mme.pca$vectors) #extract eigenvectors 
pca.vectors <- as_tibble(cbind(pop, eigenvectors)) #combine with our population assignments
df = type_convert(pca.vectors) #check all columns and convert if necessary (automatic)

#plot PC1 vs PC2 using ggplot
pca = ggplot(data = df, aes(x=V2, y=V3, fill = pop, colour = pop)) +
  geom_point(size = 5, shape = 21) +
  xlab("Principal component 1") +
  ylab("Principal component 2")

#Save plot as pdf
ggsave(filename = "/workdir/arne/physalia_lcwgs_data/data_practicals/Results/pca_pcangsd_plot.pdf", plot = pca)