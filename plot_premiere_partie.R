library(dplyr)
library(ggplot2)

rm(list=ls())
dev.off()
load("gene_expression.RData")

sum_gene_df <- as.data.frame(colSums(gene_expression)) %>% 
  rename("nb_transcrit" =`colSums(gene_expression)`) %>% 
  filter(nb_transcrit != 0)

lapoisse <- data.frame(poisson = rpois(450000, 2.18)) %>% group_by(poisson) %>% summarise(number_poisson = n())
lapoisse <- data.frame(lapoisse[2:15, 2])
lapoisse <- rbind(lapoisse, c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0))
lapoisse[is.na(lapoisse)] = 0
super <- data.frame(sum_gene_df %>% filter(nb_transcrit < 25) %>% group_by(nb_transcrit) %>% summarise(number = n()))
super <- cbind(super, poisson = lapoisse[1:24,], yend = rep(0,24))
superplot <- ggplot(data = super,
                    aes(x = nb_transcrit, y = poisson, xend = nb_transcrit, yend = yend)) +
  geom_col(stat='identity', aes(x = nb_transcrit, y = number), 
           binwidth = 2, 
           colour = "#65C6ED", fill = "#75D7FF") +
  geom_point(color = "hotpink1", size=1.5) +
  geom_segment(color = "hotpink1", linewidth=1.1) +
  theme_bw() +
  ylab("Effectif") +
  xlab("Nombre de transcrits")


ggsave(plot=superplot, filename="poisson_melange.png", width=8, height=6)



ggplot(sum_gene_df %>% filter(nb_transcrit < 25), aes(x = nb_transcrit)) +
  geom_histogram(aes(x = nb_transcrit), binwidth = 1, color = "#65C6ED", fill = "#75D7FF") +
  theme_bw()+
  ylab("Effectif") +
  xlab("Nombre de transcrits")

p <- ggplot(sum_gene_df, aes(x = seq_along(nb_transcrit), y = nb_transcrit)) +
  geom_point(shape = '.') +
  scale_y_continuous(trans = 'log10') +
  theme_bw()


p1 <- p +
  geom_hline(yintercept = 1, linetype = "dashed", color = '#75D7FF', linewidth = 0.6) +
  geom_hline(yintercept = 2000, linetype = "dashed", color = '#75D7FF', linewidth = 0.6) +
  geom_rect(aes(xmin = 0, xmax = 327392, ymin = 1, ymax = 2000),
            fill = '#75D7FF',
            alpha = 0.01,
            inherit.aes = F) +
  ylab("Nombre de transcrits") +
  xlab(" ")
p1

ggsave(plot=p1, filename="ndp2_groupe1.png", width=4, height=4)

p2 <- p +
  geom_hline(yintercept = 2000, linetype = "dashed", color = '#75D7FF', linewidth = 0.6) +
  geom_hline(yintercept = 7500, linetype = "dashed", color = '#75D7FF', linewidth = 0.6) +
  geom_rect(aes(xmin = 0, xmax = 327392, ymin = 2000, ymax = 7500),
            fill = '#75D7FF',
            alpha = 0.01,
            inherit.aes = F) +
  ylab("Nombre de transcrits") +
  xlab(" ")
ggsave(plot=p2, filename="ndp2_groupe2.png", width=4, height=4)

p3 <- p +
  geom_hline(yintercept = 5000, linetype = "dashed", color = '#75D7FF', linewidth = 0.6) +
  geom_hline(yintercept = 100000, linetype = "dashed", color = '#75D7FF', linewidth = 0.6) +
  geom_rect(aes(xmin = 0, xmax = 327392, ymin = 5000, ymax = 100000),
            fill = '#75D7FF',
            alpha = 0.01,
            inherit.aes = F) +
  ylab("Nombre de transcrits") +
  xlab(" ")
ggsave(plot=p3, filename="ndp2_groupe3.png", width=4, height=4)
