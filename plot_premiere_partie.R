library(dplyr)
library(ggplot2)
library(VaRES) #logcauchy
rm(list=ls())
dev.off()
load("gene_expression.RData")

sum_gene_df <- as.data.frame(colSums(gene_expression)) %>% 
  rename("nb_transcrit" =`colSums(gene_expression)`) %>% 
  filter(nb_transcrit != 0)

lapoisse <- data.frame(poisson = rpois(450000, 0.708)) %>% group_by(poisson) %>% summarise(number_poisson = n())
lapoisse <- rbind(lapoisse,c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0),  c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0), c(0))

labinom <- data.frame(binomial = rnbinom(77000, size = 26.9, prob = 0.73)) %>% group_by(binomial) %>% summarise(number_binomial = n())
ladeux <- data.frame(cbind(lapoisse[2:97,], labinom[5:100, ]) )

#ladeux <- data.frame(cbind(labinom[2:100, 2]))
#lapoisse[is.na(ladeux)] = 0
super <- data.frame(sum_gene_df %>% filter(nb_transcrit < 30) %>% group_by(nb_transcrit) %>% summarise(number = n()))
super <- cbind(super, poisson = ladeux[1:29,2], binomiale = ladeux[1:29,4], yend = rep(0,29))
superplot <- ggplot(data = super,
                    aes(x = nb_transcrit, xend = nb_transcrit, yend = yend)) +
  geom_col(stat='identity', aes(x = nb_transcrit, y = number), 
           binwidth = 1, 
           colour = "#65C6ED", fill = "#75D7FF") +
  geom_point(aes(y = poisson), color = "hotpink1", size=1.5) +
  geom_segment(aes(y = poisson), color = "hotpink1", linewidth=1.1) +
  geom_point(aes(y = binomiale), color = "yellow", size=1.5) +
  geom_segment(aes(y = binomiale), color = "yellow", linewidth=1.1) +
  stat_function(fun = function(x) {(dlogcauchy(x, mu = 2.18, sigma = 0.21))*150000})+
  theme_bw() +
  ylab("Effectif") +
  xlab("Nombre de transcrits")

superplot
ggsave(plot=superplot, filename="poisson_melange.png", width=8, height=6)

## On peut compter le nombre d'individu par groupe et leur proportion avant le modèle
pop_classe_avant <- count(sum_gene_df, nb_transcrit)
pop_classe_avant$proportion <- (pop_classe_avant$n / sum(pop_classe_avant$n))
test <- pop_classe_avant %>% filter(nb_transcrit < 25)
test$poisson <- c(0.69, 0.24, 0.056, 0.0098, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
test$nbinom <- dnbinom(1:10)

ggplot(test, aes(x = nb_transcrit)) +
  geom_col(aes(x = nb_transcrit, y = proportion),just = 0.5, width = 1, color = "#65C6ED", fill = "#75D7FF") +
  theme_bw()+
  geom_point(aes(y = poisson), color = "hotpink1", size=1.5) +
  ylab("Proportion") +
  xlab("Nombre de transcrits")

ggplot(pop_classe_avant %>% filter(nb_transcrit > 50)%>%filter(nb_transcrit < 50000), aes(x = nb_transcrit)) +
  geom_col(aes(x = nb_transcrit, y = proportion), color = "#65C6ED", fill = "#75D7FF") +
  theme_bw()+
  stat_function(aes(color= "distribution 1", alpha = 0.4),size = 1.1, fun = function(x) {(dnorm(x, mean = 5626, sd = 5068.9)*0.028*10)})+
  stat_function(aes(color= "distribution 2", alpha = 0.4), size = 1.1, fun = function(x) {(dnorm(x, mean = 2*5626, sd = sqrt(2)*5068.9)*0.0026*10)})+
  scale_color_manual(name = "Distributions",
                     values = c("distribution 1" = "red",
                                "distribution 2" = "darkgreen")) +
  ylab("Proportion") +
  xlab("Nombre de transcrits")

ggplot(sum_gene_df, aes(x = nb_transcrit)) +
  geom_histogram(aes(x = nb_transcrit, y = after_stat(..count../sum(..count..))), color = "#65C6ED", fill = "#75D7FF") +
  theme_bw()+
  stat_function(fun = function(x) {(dweibull(x, shape = 2.8, scale = 10.99)*0.2)})+
  ylab("Proportion") +
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
