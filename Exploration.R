library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)
data <- Read10X_h5("raw_feature_bc_matrix.D7.h5")

gene_expression <- data$`Gene Expression`
antibody_capture <- data$`Antibody Capture`

summary(gene_expression)

df <- as.data.frame(gene_expression[c(1:200), c(1:50)])

sum_gene_df <- as.data.frame(colSums(gene_expression)) %>% 
  rename("nb_transcrit" =`colSums(gene_expression)`) %>% 
  filter(nb_transcrit != 0)

ggplot(sum_gene_df, aes(x = nb_transcrit)) + 
  geom_boxplot()


boxplot(sum_gene_df)

summary(sum_gene_df)
apply(sum_gene_df)

freq_table <- as.data.frame(ftable(as.data.frame(colSums(gene_expression))))
freq_table$colSums.gene_expression.


plot(freq_table)

ggplot(freq_table, aes(x = as.numeric(colSums.gene_expression.), y = Freq)) + 
  geom_point() +
  scale_x_continuous(trans='log10') +
  geom_jitter()

