### ESSAi avec mclust

echantillon_gene <- sum_gene_df[1:200,1]
echantillon_gene <- data.frame(echantillon_gene)
library(dplyr)
plot(echantillon_gene)

library(mclust)
gene.mclust <- Mclust(echantillon_gene, G=3)
summary(gene.mclust)
plot(gene.mclust, what = "classification")
