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

fit3 <- VarSelCluster(sum_gene_df, g=3)
summary(fit3)
head(fit3@partitions@tik)
fit3@param@pi
plot(fit3, type="probs-class")


# Modèle de mélange

fit <- normalmixEM(sum_gene_df$nb_transcrit, k = 3)

classify <- function(x, means, variances, prior_probs) {
  prob <- matrix(0, nrow = 327395, ncol = 3)
  for (i in 1:327395) {
    for (j in 1:3) {
      if (j == 2) {
        prob[i, j] <- dpois(x[i], lambda = means[j]) * prior_probs[j]
      } 
      if (j == 1) {
        prob[i, j] <- dlnorm(x[i], meanlog = means[j], sdlog=variances[j]) * prior_probs[j]
      } 
      if (j == 3) {
        prob[i, j] <- dlnorm(x[i], meanlog = 2*means[j], sdlog=variances[j]) * prior_probs[j]
      } 
    }
  }
  class <- apply(prob, 1, which.max)
  return(class)
}
fit$lambda
fit$mu
fit$sigma

# Attribution des classes à chaque observation
cluster_assignment <- classify(sum_gene_df$nb_transcrit, fit$mu, fit$sigma, fit$lambda)

df_avec_classe <- data.frame(sum_gene_df, cluster_assignment)

# Modification des paramètres pour respecter les lois souhaitées
fit$lambda[1] <- fit$mu[1]
fit$sigma[2:3] <- 1

# Vérification de la qualité de l'ajustement
plot(density(sum_gene_df$nb_transcrit), xlim = range(sum_gene_df$nb_transcrit))
lines(density(dnorm(sum_gene_df$nb_transcrit, mean = fit$mu[1], sd = fit$sigma[1])), col = "red")
lines(density(dpois(sum_gene_df$nb_transcrit, lambda = fit$lambda[1])), col = "red")
lines(density(dnorm(sum_gene_df$nb_transcrit, mean = fit$mu[2], sd = fit$sigma[2])), col = "blue")
lines(density(dnorm(sum_gene_df$nb_transcrit, mean = fit$mu[3], sd = fit$sigma[3])), col = "green")