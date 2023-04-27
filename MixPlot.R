# LIBRARY ####
library(ggplot2)
library(dplyr)
library(grid)
library(tidyverse)
library(mixtools)

# DATA ####
rm(list = ls())
# Graine aléatoire
set.seed(86)
load("sum_gene_df.RData")
# Grosse triche pour le graphique (on enleve la loi de poisson des données car je sais pas faire avec)
sum_gene2 <- sum_gene_df %>%
  filter(nb_transcrit < 15000,
         nb_transcrit > 900)



# Moyenne et écartypes d'une normale qui correspondent le mieux la bosse de nos données réelles
nb_mean <- 4300
nb_sd <- 1500

# Création d'un échantillon avec les proportions correpondantes
# 2 lois normales avec mu2 = 2*mu1...
# On peut jouer sur les proportions pour faire bouger l'intensité des courbes rouge et verte
observations <- tibble(value = c(
  rnorm(n = 10000*19/28, mean = nb_mean, sd = nb_sd),
  rnorm(n = 10000*4/28, mean = 2*nb_mean, sd = sqrt(2)*nb_sd)))

# Modele de mélange pour les lois
my_mix <- normalmixEM(observations$value, k = 2)
plot(my_mix, which = 2)

# Résultats
# Moyenne
my_mix[["mu"]]
# Ecartype
my_mix[["sigma"]]
# Porportions
my_mix[["lambda"]]



# HISTOGRAMME NORMALE ####

ggplot(sum_gene2, aes(x = nb_transcrit)) +
  geom_histogram(binwidth = 200, fill = "#75D7FF", colour = "#65C6ED") +
  stat_function(
    size = 1.3,
    aes(color= "distribution 1", alpha = 0.4),
    fun = function(x) {
      (dnorm(x, mean = my_mix$mu[1], sd = my_mix$sigma[1])) *
        length(observations$value)* 200 * my_mix$lambda[1]
    }
  ) +
  stat_function(
    size = 1.3,
    aes(color= "distribution 2", alpha = 0.4),
    fun = function(x) {
      (dnorm(x, mean = my_mix$mu[2], sd = my_mix$sigma[2])) *
        length(observations$value)* 200 * my_mix$lambda[2]
    }
  ) +
  scale_color_manual(name = "Distributions",
                     values = c("distribution 1" = "red",
                                "distribution 2" = "darkgreen")) +
  ylab("Effectif") +
  xlab("Nombre de transcrits") +
  theme_bw() +
  theme(legend.position = "none")
