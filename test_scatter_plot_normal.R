library(ggplot2)

## On peut générer un échantillon de notre loi de mélange et tracer son histogramme
taille_echantillon <- 300000
echantillon_hat <- c(rpois(taille_echantillon*pi1_hat, lambda_hat), 
                     round(rnorm(taille_echantillon*pi2_hat, mu_hat, sigma_hat)), 
                     round(rnorm(taille_echantillon*pi3_hat, 2*mu_hat, sqrt(2)*sigma_hat)))

echantillon_hat <- as.data.frame(cbind(echantillon_hat, 
                                       c(rep(1, taille_echantillon*pi1_hat), 
                                         rep(2, taille_echantillon*pi2_hat), 
                                         rep(3, taille_echantillon*pi3_hat))))
colnames(echantillon_hat) <- c("echantillon", "groupe")
echantillon_hat <- echantillon_hat[which(echantillon_hat$echantillon > 0),]
echantillon_hat <- echantillon_hat[sample(nrow(echantillon_hat)),]

ggplot(echantillon_hat, aes(x = seq_along(echantillon), y = echantillon, color = as.factor(groupe))) +
  geom_point() +
  scale_y_continuous(trans = 'log10') +
  theme_bw() +
  ggtitle("Scatter plot,\n and a first determination of the groups") +
  xlab("Index") + ylab("Number of transcript") +
  scale_color_manual(values=c("black", "#66CCFF", "red"))

ggplot(sum_gene_df, aes(x = seq_along(nb_transcrit), y = nb_transcrit)) +
  geom_point(shape = '.') +
  scale_y_continuous(trans = 'log10') +
  theme_bw() +
  ggtitle("Scatter plot,\n and a first determination of the groups") +
  xlab("Index") + ylab("Number of transcript")

