################################################################
##--------------------RESULTATS OBTENUS-----------------------##
################################################################
pop_classe_avant <- count(sum_gene_df, nb_transcrit)
pop_classe_avant$proportion <- (pop_classe_avant$n / sum(pop_classe_avant$n))
proportion_25 <- pop_classe_avant %>% filter(nb_transcrit < 25)
proportion_25$poisson <- exp(dpois(1:24, 2.18, log = TRUE)  - logspace.sub(0, -2.18))

plot_poisson <- ggplot(proportion_25, aes(x = nb_transcrit)) +
  geom_col(aes(x = nb_transcrit, y = proportion),just = 0.5, width = 0.99, color = "#65C6ED", fill = "#75D7FF") +
  theme_bw()+
  geom_point(aes(y = poisson*pi1_hat, color = "Fonction de masse de la loi de Poisson"), size=1.5) +
  ylab("Proportion") +
  xlab("Nombre de transcrits") +
  scale_color_manual(values = c("Fonction de masse de la loi de Poisson" = "hotpink1"))+
  theme(legend.position = "top",
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="#909090")) +
  theme(legend.position = c(0.8, 0.94),
        legend.direction = "horizontal")+
  theme(legend.title = element_blank()) +
  labs(title = "Fonction de masse de la loi de Poisson superposées aux données réelles") +
  theme(plot.title = element_text(face = "bold",size = 13, hjust = 0, vjust = 0))


plot_poisson
ggsave(plot=plot_poisson, filename="plot_poisson.png", width=8, height=5)


sum(sum_gene_df %>% filter(nb_transcrit > 300)%>%filter(nb_transcrit < 20000))
resultats_obtenus_normales <- ggplot(sum_gene_df %>% filter(nb_transcrit > 300)%>%filter(nb_transcrit < 20000), aes(x = nb_transcrit)) +
  geom_histogram(aes(x = nb_transcrit), binwidth = 200, color = "#65C6ED", fill = "#75D7FF") +
  theme_bw()+
  stat_function(aes(color= "Première normale"),size = 1, fun = function(x) {(dnorm(x, mean = 2471, sd = 4382)*pi2_hat*90003102)})+
  stat_function(aes(color= "Seconde normale"), size = 1, fun = function(x) {(dnorm(x, mean = 2*2471, sd = sqrt(2)*4382)*pi3_hat*90003102)})+
  scale_color_manual(values = c("Première normale" = "red",
                                "Seconde normale" = "darkgreen")) +
  theme(legend.position = "top",
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="#909090")) +
  theme(legend.position = c(0.8, 0.94),
        legend.direction = "horizontal")+
  theme(legend.title = element_blank())+
  ylab("Proportion") +
  xlab("Nombre de transcrits") +
  labs(title = "Distributions des lois normales superposées aux données réelles") +
  theme(plot.title = element_text(face = "bold",size = 13, hjust = 0, vjust = 0))


resultats_obtenus_normales
ggsave(plot=resultats_obtenus_normales, filename="resultats_obtenus_normales.png", width=8, height=5)



################################################################
##--------------------RESULTATS ATTENDUS----------------------##
################################################################
pop_classe_avant <- count(sum_gene_df, nb_transcrit)
pop_classe_avant$proportion <- (pop_classe_avant$n / sum(pop_classe_avant$n))
proportion_25 <- pop_classe_avant %>% filter(nb_transcrit < 25)
proportion_25$poisson <- exp(dpois(1:24, 0.7, log = TRUE)  - logspace.sub(0, -0.7))

plot_poisson_voulu <- ggplot(proportion_25, aes(x = nb_transcrit)) +
  geom_col(aes(x = nb_transcrit, y = proportion),just = 0.5, width = 0.99, color = "#65C6ED", fill = "#75D7FF") +
  theme_bw()+
  geom_point(aes(y = poisson*pi1_hat, color = "Fonction de masse de la loi de Poisson"), size=1.5) +
  ylab("Proportion") +
  xlab("Nombre de transcrits") +
  scale_color_manual(values = c("Fonction de masse de la loi de Poisson" = "hotpink1"))+
  theme(legend.position = "top",
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="#909090")) +
  theme(legend.position = c(0.8, 0.94),
        legend.direction = "horizontal")+
  theme(legend.title = element_blank())+
  labs(title = "Fonction de masse d'une loi de Poisson superposées aux données réelles") +
  theme(plot.title = element_text(face = "bold",size = 13, hjust = 0, vjust = 0))
  
plot_poisson_voulu
ggsave(plot=plot_poisson_voulu, filename="plot_poisson_voulu.png", width=8, height=5)


normale_voulu <- ggplot(sum_gene_df %>% filter(nb_transcrit > 300)%>%filter(nb_transcrit < 20000), aes(x = nb_transcrit)) +
  geom_histogram(aes(x = nb_transcrit), binwidth = 200, color = "#65C6ED", fill = "#75D7FF") +
  theme_bw()+
  stat_function(aes(color= "Première normale"),size = 1, fun = function(x) {(dnorm(x, mean = 4400, sd = 1600)*pi2_hat*51033102)})+
  stat_function(aes(color= "Seconde normale"), size = 1, fun = function(x) {(dnorm(x, mean = 2*4400, sd = sqrt(2)*1600)*pi3_hat*100933102)})+
  scale_color_manual(values = c("Première normale" = "red",
                                "Seconde normale" = "darkgreen")) +
  theme(legend.position = "top",
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="#909090")) +
  theme(legend.position = c(0.8, 0.94),
        legend.direction = "horizontal")+
  theme(legend.title = element_blank())+
  ylab("Proportion") +
  xlab("Nombre de transcrits") +
  labs(title = "Distributions de deux lois normales superposées aux données réelles") +
  theme(plot.title = element_text(face = "bold",size = 13, hjust = 0, vjust = 0))


normale_voulu
ggsave(plot=normale_voulu, filename="normale_voulu.png", width=8, height=5)


################################################################
##--------------------RESULTATS WEIBULL----------------------##
################################################################
pop_classe_avant <- count(sum_gene_df, nb_transcrit)
pop_classe_avant$proportion <- (pop_classe_avant$n / sum(pop_classe_avant$n))
proportion_25 <- pop_classe_avant %>% filter(nb_transcrit < 25)
proportion_25$poisson <- exp(dpois(1:24, lambda_hat, log = TRUE)  - logspace.sub(0, -lambda_hat))

plot_poisson_weibull <- ggplot(proportion_25, aes(x = nb_transcrit)) +
  geom_col(aes(x = nb_transcrit, y = proportion),just = 0.5, width = 0.99, color = "#65C6ED", fill = "#75D7FF") +
  theme_bw()+
  geom_point(aes(y = poisson*pi1_hat, color = "Fonction de masse de la loi de Poisson"), size=1.5) +
  stat_function(aes(color= "Densité de la loi de Weibull"),size = 1, fun = function(x) {(dweibull(x, shape = alpha_hat, scale = beta_hat)*pi4_hat)})+
  ylab("Proportion") +
  xlab("Nombre de transcrits") +
  scale_color_manual(values = c("Fonction de masse de la loi de Poisson" = "hotpink1", "Densité de la loi de Weibull" = "#efd970"))+
  theme(legend.position = "top",
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="#909090")) +
  theme(legend.position = c(0.65, 0.94),
        legend.direction = "horizontal")+
  theme(legend.title = element_blank())+
  labs(title = "Fonction de masse d'une loi de Poisson superposées aux données réelles") +
  theme(plot.title = element_text(face = "bold",size = 13, hjust = 0, vjust = 0))

plot_poisson_weibull
ggsave(plot=plot_poisson_weibull, filename="plot_poisson_weibull.png", width=8, height=5)


normale_weibull <- ggplot(sum_gene_df %>% filter(nb_transcrit > 300)%>%filter(nb_transcrit < 20000), aes(x = nb_transcrit)) +
  geom_histogram(aes(x = nb_transcrit), binwidth = 200, color = "#65C6ED", fill = "#75D7FF") +
  theme_bw()+
  stat_function(aes(color= "Première normale"),size = 1, fun = function(x) {(dnorm(x, mean = mu_hat, sd = sigma_hat)*pi2_hat*161033102)})+
  stat_function(aes(color= "Seconde normale"), size = 1, fun = function(x) {(dnorm(x, mean = 2*mu_hat, sd = sqrt(2)*sigma_hat)*pi3_hat*161033102)})+
  scale_color_manual(values = c("Première normale" = "red",
                                "Seconde normale" = "darkgreen")) +
  theme(legend.position = "top",
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="#909090")) +
  theme(legend.position = c(0.8, 0.94),
        legend.direction = "horizontal")+
  theme(legend.title = element_blank())+
  ylab("Proportion") +
  xlab("Nombre de transcrits") +
  labs(title = "Distributions des deux lois normales superposées aux données réelles") +
  theme(plot.title = element_text(face = "bold",size = 13, hjust = 0, vjust = 0))


normale_weibull
ggsave(plot=normale_weibull, filename="normale_weibull.png", width=8, height=5)
