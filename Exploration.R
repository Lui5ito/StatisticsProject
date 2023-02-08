library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)

#On récupère les données fournies
data <- Read10X_h5("raw_feature_bc_matrix.D7.h5")

#Les données sont sous formes d'une "liste" de matrices. On les extraits. 
gene_expression <- data$`Gene Expression`
antibody_capture <- data$`Antibody Capture`

#Mais ce sont encore des objet particulier. On ne peut pas les trasnformer en dataframe car les matrices sont trop grosses.
summary(gene_expression)
str(gene_expression)
summary(gene_expression@i)
summary(gene_expression@p)
summary(gene_expression@Dim)
summary(gene_expression@x)
summary(gene_expression@Dimnames)
summary(gene_expression@factors)

#####On peut commencer par extraire une 'sous' matrice de gene_expression pour l'analyser. Mais il va falloir trouver une manière de faire des GLM, modèles de mélanges sur des dgCMatrix. Voir package 'Matrix' et 'glmnet'.
#Pour l'instant la matrice qui nous intéresse c'est gene_expression. C'est la qu'on trouve l'information des goutelettes.

gene_expression_df <- as.data.frame(gene_expression[c(1:200), c(1:50)])

summary(gene_expression_df)

#Ici on veut calculer le nombre de transcrit par goutelette.
sum_gene_df <- as.data.frame(colSums(gene_expression)) %>% 
  rename("nb_transcrit" =`colSums(gene_expression)`) %>% 
  filter(nb_transcrit != 0)

summary(sum_gene_df)

#Marche pas
ggplot(sum_gene_df, aes(x = nb_transcrit)) + 
  geom_boxplot()
boxplot(sum_gene_df)
apply(sum_gene_df)
#aMarche pas

#On veut récuperer les fréquences des goutelettes qui ont x transcrits
freq_table <- as.data.frame(ftable(as.data.frame(colSums(gene_expression))))

plot(freq_table)

#On plot sur une échelle logarithmique le nombre de goutelettes qui ont le même nombre de transcrit.
#Logarithmique car le nombre de transcrit commence à 0 et peut aller à plusieurs milliers.
ggplot(freq_table, aes(x = as.numeric(colSums.gene_expression.), y = Freq)) + 
  geom_point() +
  scale_x_continuous(trans='log10') +
  geom_jitter()

