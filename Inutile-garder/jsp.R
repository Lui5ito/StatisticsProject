library(dplyr)
library(ggplot2)

rm(list=ls())
dev.off()
load("sum_gene_df.RData")

lapoisse <- data.frame(poisson = rpois(450000, 2.18)) %>% group_by(poisson) %>% summarise(number_poisson = n())
lapoisse <- data.frame(lapoisse[2:15, 2])
lapoisse <- rbind(lapoisse, c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0))
lapoisse[is.na(lapoisse)] = 0
super <- data.frame(sum_gene_df %>% filter(nb_transcrit < 25) %>% group_by(nb_transcrit) %>% summarise(number = n()))
super <- cbind(super, poisson = lapoisse[1:24,], yend = rep(0,24))
superplot <- ggplot(data = super,
                    aes(x = nb_transcrit, y = poisson, xend = nb_transcrit, yend = yend)) +
  geom_col(stat='identity', aes(x = nb_transcrit, y = number, color = "Données réelles"),
           binwidth = 2,
           colour = "#65C6ED", fill = "#75D7FF") +
  geom_point(color = "hotpink1", size=1.5) +
  geom_segment(linewidth=1, aes(color = "Simulation groupe 1")) +
  theme_bw() +
  ylab("Effectif") +
  xlab("Nombre de transcrits") +
  scale_color_manual(name = "Distributions",
                     values = c("Données réelles" = "#75D7FF",
                                "Simulation groupe 1" = "hotpink1")) +
  theme(legend.position = "top",
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="#909090")) +
  theme(legend.title = element_blank())


superplot
ggsave(plot=superplot, filename="poisson_melange.png", width=8, height=6)
