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
library(mixtools)
fit <- normalmixEM(sum_gene_df$nb_transcrit, k = 3)

classify_avec_moyenne_du_modèle <- function(x, means, variances, prior_probs) {
  prob <- matrix(0, nrow = 327395, ncol = 3)
  for (i in 1:327395) {
    #On calcul la proba d'être dans la classe 2, ie. de suivre une loi de Poisson.
    prob[i, 2] <- dpois(x[i], lambda = means[2]) * prior_probs[2]
    
    #On calcul la proba d'être dans la classe 1, ie. de suivre une normale de moyenne n
    prob[i, 1] <- dlnorm(x[i], meanlog = means[1], sdlog=variances[1]) * prior_probs[1]
    
    #On calcul la proba d'être dans la classe 3, ie. de suivre une normale de moyenne 2n idéalement.
    prob[i, 3] <- dlnorm(x[i], meanlog = 2*means[3], sdlog=variances[3]) * prior_probs[3]
    }
  class <- apply(prob, 1, which.max)
  return(class)
}

classify_avec_moyenne_voulue_2n <- function(x, means, variances, prior_probs) {
  prob <- matrix(0, nrow = 327395, ncol = 3)
  for (i in 1:327395) {
    #On calcul la proba d'être dans la classe 2, ie. de suivre une loi de Poisson.
    prob[i, 2] <- dpois(x[i], lambda = means[2]) * prior_probs[2]
    
    #On calcul la proba d'être dans la classe 1, ie. de suivre une normale de moyenne n
    prob[i, 1] <- dlnorm(x[i], meanlog = means[1], sdlog=variances[1]) * prior_probs[1]
    
    #On calcul la proba d'être dans la classe 3, ie. de suivre une normale de moyenne 2n idéalement.
    prob[i, 3] <- dlnorm(x[i], meanlog = 2*means[1], sdlog=variances[3]) * prior_probs[3]
  }
  class <- apply(prob, 1, which.max)
  return(class)
}
fit$lambda
fit$mu
fit$sigma

# Attribution des classes à chaque observation
cluster_assignment_avec_moyenne_du_modèle <- classify_avec_moyenne_du_modèle(sum_gene_df$nb_transcrit, fit$mu, fit$sigma, fit$lambda)
df_avec_classe_avec_moyenne_du_modèle <- data.frame(sum_gene_df, cluster_assignment_avec_moyenne_du_modèle)

summary(df_avec_classe_avec_moyenne_du_modèle)
ggplot(df_avec_classe_avec_moyenne_du_modèle, aes(x = as.factor(cluster_assignment_avec_moyenne_du_modèle), y = nb_transcrit)) +
  geom_violin(aes(fill = as.factor(cluster_assignment_avec_moyenne_du_modèle)), trim = FALSE) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"))


cluster_assignment_avec_moyenne_voulue_2n <- classify_avec_moyenne_voulue_2n(sum_gene_df$nb_transcrit, fit$mu, fit$sigma, fit$lambda)
df_avec_classe_avec_moyenne_voulue_2n <- data.frame(sum_gene_df, cluster_assignment_avec_moyenne_voulue_2n)

summary(df_avec_classe_avec_moyenne_voulue_2n)
ggplot(df_avec_classe_avec_moyenne_voulue_2n, aes(x = as.factor(cluster_assignment_avec_moyenne_voulue_2n), y = nb_transcrit)) +
  geom_violin(aes(fill = as.factor(cluster_assignment_avec_moyenne_voulue_2n)), trim = FALSE) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"))


range(df_avec_classe$nb_transcrit)

# Vérification de la qualité de l'ajustement
plot(density(sum_gene_df$nb_transcrit), xlim = range(sum_gene_df$nb_transcrit))
lines(density(dnorm(sum_gene_df$nb_transcrit, mean = fit$mu[1], sd = fit$sigma[1])), col = "red")
lines(density(dpois(sum_gene_df$nb_transcrit, lambda = fit$lambda[1])), col = "red")
lines(density(dnorm(sum_gene_df$nb_transcrit, mean = fit$mu[2], sd = fit$sigma[2])), col = "blue")
lines(density(dnorm(sum_gene_df$nb_transcrit, mean = fit$mu[3], sd = fit$sigma[3])), col = "green")