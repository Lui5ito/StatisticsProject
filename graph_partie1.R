# LIBRARY ####
library(ggplot2)
library(dplyr)
library(grid)
library(plotly)

# HISTOGRAMME NORMALE ####
ggplot(sum_gene_df, aes(x = nb_transcrit)) +
  geom_histogram(aes(y = ..density..), fill = "#66CCFF", color = "blue", bins = 100) +
  xlim(1200, 9000) +
  ylim(0, 0.0004) +
  theme_bw() +
  stat_function(fun = dnorm, args = list(mean = 4350, sd = 1450), color = "red", size = 1)


ggplotly(gg)

hist(sum_gene_df$nb_transcrit, freq = FALSE, col = "#66CCFF", ylim = c(0, 0.0004), xlim = c(1200, 9000))

# Ajout de la courbe de densité empirique
dens <- density(sum_gene_df$nb_transcrit)
lines(dens, col = "red", lwd = 2)
