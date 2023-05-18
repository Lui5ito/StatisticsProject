setwd("~/2A/statproj1")
# LIBRARY ####
library(ggplot2)
library(dplyr)
library(grid)
library(plotly)
library(tidyverse)
library(mixtools)

# DATA ####
rm(list = ls())
# Graine aléatoire
set.seed(86)
load("sum_gene_df.RData")
load("gene_expression.RData")

sum_gene_df <- as.data.frame(colSums(gene_expression)) %>% 
  rename("nb_transcrit" =`colSums(gene_expression)`) %>% 
  filter(nb_transcrit != 0)

group_by(nb_transcrit)%>%
  summarize(n=n())%>%
  mutate(prop=n/sum(n))

# Grosse triche pour le graphique (on enleve la loi de poisson des données car je sais pas faire avec)
sum_gene2 <- sum_gene_df %>%
  filter(nb_transcrit < 40000,
         nb_transcrit > 200)




ggplot(sum_gene2, aes(x = nb_transcrit)) +
  geom_histogram(binwidth = 20,fill = "#75D7FF", colour = "#65C6ED", aes(color = "Données réelles"))

plot(1:100, dnbinom(1:100, size = 10, prob = 0.8))  


# Moyenne et écartypes d'une normale qui correspondent le mieux la bosse de nos données réelles
nb_mean <- 5787
nb_sd <- 4913.9

# Création d'un échantillon avec les proportions correpondantes
# 2 lois normales avec mu2 = 2*mu1...
# On peut jouer sur les proportions pour faire bouger l'intensité des courbes rouge et verte
observations <- tibble(value = c(
  rnorm(n = 10000*(0.028)/(0.028+0.003), mean = nb_mean, sd = nb_sd),
  rnorm(n = 10000*(0.003)/(0.028+0.003), mean = 2*nb_mean, sd = sqrt(2)*nb_sd)))

# Modele de mélange pour les lois
my_mix <- normalmixEM(observations$value, k = 2)
plot(my_mix, which = 3)

# Résultats
# Moyenne
my_mix[["mu"]]
# Ecartype
my_mix[["sigma"]]
# Porportions
my_mix[["lambda"]]



# HISTOGRAMME NORMALE ####

Francky <- ggplot(sum_gene2, aes(x = nb_transcrit)) +
  geom_histogram(binwidth = 100, fill = "#75D7FF", colour = "#65C6ED", aes(color = "Données réelles"))+
  stat_function(
    size = 1.1,
    aes(color= "Densité du groupe 2"),
    fun = function(x) {
      (dnorm(x, mean = my_mix$mu[1], sd = my_mix$sigma[1])) *
        length(observations$value)* 200 * my_mix$lambda[1]
    }
  ) +
  stat_function(
    size = 1.1,
    aes(color= "Densité du groupe 3"),
    fun = function(x) {
      (dnorm(x, mean = my_mix$mu[2], sd = my_mix$sigma[2])) *
        length(observations$value)* 200 * my_mix$lambda[2]
    }
  ) +
  scale_color_manual(values = c("Données réelles" = "#75D7FF",
                                "Densité du groupe 2" = "red",
                                "Densité du groupe 3" = "darkgreen"))+
  ylab("Densité") +
  xlab("Nombre de transcrits") +
  theme_bw() +
  theme(legend.position = "top",
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="#909090")) +
  theme(legend.title = element_blank())

Francky
ggsave(plot = Francky,
       filename = "hist_froc.png",
       width = 8,
       height = 6)
Francky