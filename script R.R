rm(list=ls())

### PACKAGES NECESSAIRES ###
install.packages("Seurat")
install.packages("hdf5r")
install.packages("rhdf5")
install.packages("mclust")
install.packages("VarSelLCM")
library(Seurat)
require(VarSelLCM)

### IMPORT FICHIER ### 
brec <- Read10X_h5("C:/Users/julie/Documents/ENSAI/Année 2/s2/Projet stat/BREC/raw_feature_bc_matrix.D7.h5", use.names = TRUE, unique.features = TRUE)
gene_expression <- brec[[1]]
antibody_capture <- brec[[2]]
colnames(antibody_capture)
dimnames(antibody_capture)
dimnames(gene_expression)
gene_expression[1:10,1]

df = as.data.frame(gene_expression[,1:100])
df_sum <- colSums(df)
nom_goutellettes <- as.data.frame(gene_expression[0,1:100])
nouveau_df <- rbind(nom_goutellettes, df_sum)





install.packages("mixtools")
library(mixtools)

# Construire un modèle de mélange à trois composantes
# Composante 1: loi de Poisson pour les gouttelettes vides
# Composante 2: loi log-normale pour les gouttelettes contenant une cellule
# Composante 3: loi log-normale d'espérance 2n pour les gouttelettes contenant deux cellules
model <- normalmixEM(df_sum, k = 3, lambda = c(1, 1, 1),
                     mu = c(0, log(mean(df_sum)), log(2*mean(df_sum))),
                     sigma = c(1, 1, 1))

# Afficher les paramètres estimés pour chaque composante
print(model$mu)
print(model$sigma)
print(model$lambda)
#mu est la moyenne de la distribution, c'est la valeur autour de laquelle les données sont censées être regroupées.
#sigma est l'écart-type de la distribution, il mesure la dispersion des données autour de la moyenne.
#lambda est la proportion de chaque composante dans le mélange. Il mesure la fréquence de chaque composante dans les données

# Prévoir les classes pour chaque gouttelette
class_predictions <- predict(model, fdef, df_sum)


# Installer et charger la librairie mclust
install.packages("mclust")
library(mclust)

# Construire le modèle de mélange
fit2 <- Mclust(df_sum, G = 3) # On spécifie le nombre de composantes à 3 ici

# Afficher les résultats
summary(fit2)
plot(fit2)
plot(fit)

# Prédire la classe pour chaque observation
classe <- predict(fit)

# Afficher les probabilités latentes pour chaque observation
probas <- predict(fit, what = "probabilities")
