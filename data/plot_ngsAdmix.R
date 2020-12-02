library(tidyverse) #load the tidyverse package for formatting and plotting

#Load the covariance matrix
cov = read_table("/workdir/arne/physalia_lcwgs_data/data_practicals/Results/MME_ngsAdmix_K2_out.qopt", col_names = F)

#We will also add a column with population assingments
pop <- c("JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA"
         ,"PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY"
         ,"MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS"
         ,"MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU")

cov.id = as.data.frame(cbind(pop, cov))
names(cov.id) = c("pop","q1","q2")
barplot(t(as.matrix(subset(cov.id, select=q1:q2))), col=c("firebrick","royalblue"), border=NA)

#Save plot as pdf
ggsave(filename = "/workdir/arne/physalia_lcwgs_data/data_practicals/Results/pca_pcangsd_plot.pdf", plot = pca)
