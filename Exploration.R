library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)

#On récupère les données fournies
#data <- Read10X_h5("raw_feature_bc_matrix.D7.h5")

#Les données sont sous formes d'une "liste" de matrices. On les extraits. 
#gene_expression <- data$`Gene Expression`
#antibody_capture <- data$`Antibody Capture`
#save(gene_expression, file = "gene_expression.RData")
#save(antibody_capture, file = "antibody_capture.RData")

load("gene_expression.RData")
load("antibody_capture.RData")

#Mais ce sont encore des objet particulier. On ne peut pas les trasnformer en dataframe car les matrices sont trop grosses.
#Ici ça explique comment gérer un objetdgCMatrix : https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/

as.data.frame(gene_expression[c(1:3), c(1:3)])
summary(gene_expression)
str(gene_expression)

#Liste de toutes les non-zéros valeures.
gene_expression@x
summary(gene_expression@x)

#Renvoie le numéro de la ligne de chaque valeure non-zéro. (Correspondance avec x)
gene_expression@i
summary(gene_expression@i)

#Liste de taille le nombre de colonne, et renvoie le nombre de non-zéro par colonne (en faisant la différence, cf le site)
gene_expression@p
summary(gene_expression@p)

#Nombre de lignes et nombre de colonnes
gene_expression@Dim

#Donne le nom des lignes s'ils existent et le nom des colonnes s'ils existent
gene_expression@Dimnames

# A quoi ressemble les lignes de genes_expression
as.data.frame(gene_expression@Dimnames[[1]]) %>% slice(1:50) %>% View()
nrow_genes <- gene_expression@Dimnames[[1]] %>% length() %>% as.numeric()
# A quoi ressemble les colonnes de genes_expression
as.data.frame(gene_expression@Dimnames[[2]]) %>% slice(1:50) %>% View()
ncol_genes <- gene_expression@Dimnames[[2]] %>% length() %>% as.numeric()


#Pour l'instant la matrice qui nous intéresse c'est gene_expression. C'est la qu'on trouve l'information des goutelettes.


#Ici on veut calculer le nombre de transcrit par goutelette.
sum_gene_df <- as.data.frame(colSums(gene_expression)) %>% 
  rename("nb_transcrit" =`colSums(gene_expression)`) %>% 
  filter(nb_transcrit != 0)

summary(sum_gene_df)

#On veut récuperer les fréquences des goutelettes qui ont x transcrits
freq_table <- as.data.frame(ftable(as.data.frame(colSums(gene_expression))))
colnames(freq_table)[colnames(freq_table) == 'colSums.gene_expression.'] <- 'Nombre_de_goutelettes_qui_ont_le_même_nombre_de_transcrits'
colnames(freq_table)[colnames(freq_table) == 'Freq'] <- 'Nombre_de_transcrits'

#On plot sur une échelle logarithmique le nombre de goutelettes qui ont le même nombre de transcrit.
#Logarithmique car le nombre de transcrit commence à 0 et peut aller à plusieurs milliers.
ggplot(freq_table, aes(x = as.numeric(Nombre_de_goutelettes_qui_ont_le_même_nombre_de_transcrits), y = Nombre_de_transcrits)) + 
  geom_point() +
  scale_x_continuous(trans='log10') +
  geom_jitter() +
  labs(title = 'Nombre de transcrits par goutelettes en fonctions des goutelettes qui ont le même nombre de transcrits',subtitle = 'Echelle logarithmique en abscisse', x = 'Nombre de goutelettes qui ont le même nombre de transcrit.', y = 'Nombre de transcrits dans ces goutelettes')

