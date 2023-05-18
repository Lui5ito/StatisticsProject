library(Seurat)
library(hdf5r)

## On récupère les données fournies.
data <- Read10X_h5("raw_feature_bc_matrix.D7.h5")

## Les données sont sous formes d'une "liste" de matrices. On les extraits. 
gene_expression <- data$`Gene Expression`
antibody_capture <- data$`Antibody Capture`

## On sauvegarde les données.

save(gene_expression, file = "gene_expression.RData")
save(antibody_capture, file = "antibody_capture.RData")