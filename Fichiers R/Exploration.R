library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)

################################################################################
##------------------------------INSTRUCTIONS----------------------------------##
################################################################################
# Ce fichier contient l'import des données brutes et la sauvegarde par morceaux de celles-ci.

# Il y a aussi les premiers nuages de points.

################################################################################
##-------------------------------EXPLORATION----------------------------------##
################################################################################

## On récupère les données fournies.
data <- Read10X_h5("raw_feature_bc_matrix.D7.h5")

## Les données sont sous formes d'une "liste" de matrices. On les extraits. 
gene_expression <- data$`Gene Expression`
antibody_capture <- data$`Antibody Capture`

## On sauvegarde les données.
save(gene_expression, file = "gene_expression.RData")
save(antibody_capture, file = "antibody_capture.RData")

## On ne considère pas les gouttelettes à 0 transcrits, et on somme le nombre de transcrits par gouttelettes.

sum_gene_df <- as.data.frame(colSums(gene_expression)) %>% 
  rename("nb_transcrit" =`colSums(gene_expression)`) %>% 
  filter(nb_transcrit != 0)

## Nuage de points des données
nuage_points_donnees <- ggplot(sum_gene_df, aes(x = seq_along(nb_transcrit), y = nb_transcrit)) +
  geom_point(shape = '.') +
  scale_y_continuous(trans = 'log10') +
  theme_bw() +
  ggtitle("Nuage de points des données") +
  xlab("Index de la gouttelettes") + ylab("Nombre de transcrits")

nuage_points_donnees
ggsave(plot=nuage_points_donnees, filename="images/nuage_points_donnees.png", width=6, height=6)

