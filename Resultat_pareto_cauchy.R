library(dplyr) #nécessaire pour l'importation des données
library(DPQ) #logspace.sub
library(nloptr) #nolptr
library(doParallel) #créer les clusters
library(foreach) #permet le parallélisme
library(tictoc) #calcule le temps
library(ggplot2) #plots
library(distributionsrd) #pareto

rm(list=ls())
set.seed(1664)
################################################################################
##--------------------------------DONNEES-------------------------------------##
################################################################################
load("gene_expression.RData")

sum_gene_df <- as.data.frame(colSums(gene_expression)) %>% 
  rename("nb_transcrit" =`colSums(gene_expression)`) %>% 
  filter(nb_transcrit != 0)

sum_gene_df_sup_1000 <- sum_gene_df %>% filter(nb_transcrit > 1000)


data_sup_1000 <- sum_gene_df_sup_1000[,1]

data <- sum_gene_df[,1]
################################################################################
##-------------------------------FONCTIONS------------------------------------##
################################################################################

logvraissemblance <- function(params){
  logphi1 <- dpareto(data, k = params[3], xmin = 1, log = TRUE)
  
  logphi2 <- logspace.sub(pcauchy(data+1/2, location = params[1], scale = params[2], log.p = TRUE), pcauchy(data-1/2, location = params[1], scale = params[2], log.p = TRUE))
  logphi3 <- logspace.sub(pcauchy(data+1/2, location = 2*params[1], scale = sqrt(2)*params[2], log.p = TRUE), pcauchy(data-1/2, location = 2*params[1], scale = sqrt(2)*params[2], log.p = TRUE))
  
  
  ln1 <- log(pi1_r) + logphi1
  ln2 <- log(pi2_r) + logphi2
  ln3 <- log(pi3_r) + logphi3
  
  somme_phi_pondere <- pi1_r*exp(logphi1) + pi2_r*exp(logphi2) + pi3_r*exp(logphi3)
  
  t1 <- pi1_r*exp(logphi1) / somme_phi_pondere
  t2 <- pi2_r*exp(logphi2) / somme_phi_pondere
  t3 <- pi3_r*exp(logphi3) / somme_phi_pondere
  
  lv1 <- sum(t1*ln1)
  lv2 <- sum(t2*ln2)
  lv3 <- sum(t3*ln3)
  
  return(lv1+lv2+lv3)
}

logvraissemblance_pour_nloptr <- function(params){return(-logvraissemblance(params))}


################################################################################
##-------------------------------ALGORITHME-----------------------------------##
################################################################################

## Initialisation du parallélisme
nbre_coeurs_disponibles = detectCores()
pourcentage_des_coeurs_voulu <- 30
nbre_coeurs_voulu <- makeCluster(pourcentage_des_coeurs_voulu)
registerDoParallel(nbre_coeurs_voulu)

nbre_random_start <- 30 #### On prend plus que 2x le nombre de paramètres car on a un échantillon dispoportionné
nbre_parametre_interet <- 7 #### log-vraisemblance complétée, mu, sigma, pi1, pi2, pi3, alpha


## Initialisation de tous ce qui est statique
n <- length(data)
n_sup_1000 <- length(data_sup_1000)

tic("Application")
## Cette boucle random start EST l'algorithme EM
## resultats est une liste contenant chaque random start. 
## Un random start est une liste contenant des vecteurs qui contiennent eux même les itérations de chaque paramètre d'interêt.
resultats <- foreach (i=1:nbre_random_start, .packages=c("DPQ", "nloptr", "distributionsrd")) %dopar% {
  
  ## On setseed à chaque tour pour changer les initilisations.
  set.seed(i)
  
  ## On initialise la liste qui va contenir le résultat du random start de l'algorithme 
  ma_liste <- NULL
  
  ## (ré-)Initialisation des suivis des paramètres d'interets
  suite_mu <- c()
  suite_sigma <- c()
  suite_pi1 <- c()
  suite_pi3 <- c()
  suite_pi2 <- c()
  suite_alpha <- c()
  Lvc <- c()
  
  
  ## Initialisation aléatoire des paramètres
  mu_r <- data_sup_1000[sample(1:n_sup_1000, 1)] 
  sigma_r <- data_sup_1000[sample(1:n_sup_1000, 1)]
  alpha_r <- data[sample(1:n, 1)]
  pi1_r <- 1/3
  pi2_r <- 1/3
  pi3_r <- 1/3
  
  ## Première instance des paramètres d'interets
  suite_sigma <- c(sigma_r)
  suite_mu <- c(mu_r)
  suite_pi1 <- c(pi1_r)
  suite_pi2 <- c(pi2_r)
  suite_pi3 <- c(pi3_r)
  suite_alpha <- c(alpha_r)
  
  t <- try(repeat{
    
    ## On a bseoin de ces calculs pour lambda_r
    logphi1 <- dpareto(data, k = alpha_r, xmin = 1, log = TRUE)
    logphi2 <- logspace.sub(pcauchy(data+1/2, location = mu_r, scale = sigma_r, log.p = TRUE), pcauchy(data-1/2, location = mu_r, scale = sigma_r, log.p = TRUE))
    logphi3 <- logspace.sub(pcauchy(data+1/2, location = 2*mu_r, scale = sqrt(2)*sigma_r, log.p = TRUE), pcauchy(data-1/2, location = 2*mu_r, scale = sqrt(2)*sigma_r, log.p = TRUE))
    
    ln1 <- log(pi1_r) + logphi1
    ln2 <- log(pi2_r) + logphi2
    ln3 <- log(pi3_r) + logphi3
    
    somme_phi_pondere <- pi1_r*exp(logphi1) + pi2_r*exp(logphi2) + pi3_r*exp(logphi3)
    
    t1 <- pi1_r*exp(logphi1) / somme_phi_pondere
    t2 <- pi2_r*exp(logphi2) / somme_phi_pondere
    t3 <- pi3_r*exp(logphi3) / somme_phi_pondere
    
    lv1 <- sum(t1*ln1)
    lv2 <- sum(t2*ln2)
    lv3 <- sum(t3*ln3)
    
    T1 <- sum(t1)
    T2 <- sum(t2)
    T3 <- sum(t3)
    
    ## On doit calculer la log-vraisemblance à chaque itération puisque c'est notre porte de sortie de la boucle repeat
    LV_r <- lv1+lv2+lv3
    ## Sauvegarde de la valeur de la log-vraisemblance complétée à cette instant
    Lvc <- c(Lvc, LV_r)
    
    ## On optimise avec l'algorithme Bobyqat
    res_nlopt <- nloptr(x0 = c(mu_r, sigma_r, alpha_r), eval_f = logvraissemblance_pour_nloptr, opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 10000, tol_rel=1e-15, xtol_abs=1e-15), lb = c(1, 1, 0))
    
    ## Mises à jours des paramètres
    mu_r <- res_nlopt$solution[1]
    sigma_r <- res_nlopt$solution[2]
    alpha_r <- res_nlopt$solution[3]
    pi1_r <- T1/n
    pi2_r <- T2/n
    pi3_r <- T3/n
    
    ##Sauvegarde des historiques des paramètres d'interes
    suite_mu <- c(suite_mu, mu_r)
    suite_sigma <- c(suite_sigma, sigma_r)
    suite_pi1 <- c(suite_pi1, pi1_r)
    suite_pi2 <- c(suite_pi2, pi2_r)
    suite_pi3 <- c(suite_pi3, pi3_r)
    suite_alpha <- c(suite_alpha, alpha_r)
    
    
    
    if (length(Lvc)>2 & abs((tail(Lvc, 1) - tail(Lvc, 2)[1])) < 0.1) { break }
  })
  
  
  ## On vérifie qu'il n'y a pas eu une seule erreur de nlopt dans le repeat
  if(!(inherits(t, "try-error"))){
    
    ## On sauvegarde alors tous les paramètres d'intérêts finaux de notre algorithme EM
    Lvc <- c(Lvc, logvraissemblance(c(mu_r, sigma_r, alpha_r))) #Pour avoir autant de points que les paramètres
    ma_liste <- list(logvraisemblance  = Lvc, mu = suite_mu, sigma = suite_sigma, pi1 = suite_pi1, pi2 = suite_pi2, pi3 = suite_pi3, alpha = suite_alpha)
  }
  return (ma_liste)
}
toc()
stopCluster(nbre_coeurs_voulu)
summary(resultats)
saveRDS(resultats, file = "resultats_application_pareto_cauchy.RData")
################################################################################
##--------------------------------RESULTATS-----------------------------------##
################################################################################

resultats <- readRDS("resultats_application_pareto_cauchy.RData")

## On veut maintenant récupérer le meilleur des random starts, celui qui a la log-vraisemblance complétée la plus élevée.
max_index <- 1
length(resultats[[8]])
max_lvc <- tail(resultats[[1]][[1]], 1)
for (i in 2:length(resultats)){
  if (length(resultats[[i]]) != 0){
    if (tail(resultats[[i]][[1]], 1)>max_lvc) {
      max_lvc <- tail(resultats[[i]][[1]], 1)
      max_index <- i
    }
  }
}

################################################################
##--------------Evolution de lvc et paramètres----------------##
################################################################

## log-vraisemblance complétée
resultats[[max_index]][[1]] 
vrai <- ggplot(as.data.frame(resultats[[max_index]][[1]][1:16]), aes(x=seq_along(resultats[[max_index]][[1]][1:16]), y=resultats[[max_index]][[1]][1:16])) +
  geom_point(color = "#21ADE5", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab("log-vraisemblance complétée")
ggsave(plot=vrai, filename="iteration_lvc.png", width=8, height=5)

## mu
resultats[[max_index]][[2]]
mu <- ggplot(as.data.frame(resultats[[max_index]][[2]]), aes(x=seq_along(resultats[[max_index]][[2]]), y=resultats[[max_index]][[2]])) +
  geom_point(color = "#21ADE5", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab(expression(mu))
ggsave(plot=mu, filename="iteration_mu.png", width=8, height=5)

## sigma
resultats[[max_index]][[3]] ## pour sigma
sig <- ggplot(as.data.frame(resultats[[max_index]][[3]]), aes(x=seq_along(resultats[[max_index]][[3]]), y=resultats[[max_index]][[3]])) +
  geom_point(color = "#21ADE5", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab(expression(sigma))
ggsave(plot=sig, filename="iteration_sigma.png", width=8, height=5)

## alpha
resultats[[max_index]][[7]]
alpha <- ggplot(as.data.frame(resultats[[max_index]][[7]]), aes(x=seq_along(resultats[[max_index]][[7]]), y=resultats[[max_index]][[7]])) +
  geom_point(color = "#21ADE5", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab(expression(alpha))


## pi1
resultats[[max_index]][[4]]
ggplot(as.data.frame(resultats[[max_index]][[4]]), aes(x=seq_along(resultats[[max_index]][[4]]), y=resultats[[max_index]][[4]])) +
  geom_point(color = "#21ADE5", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab(expression(pi[1]))

## pi2
resultats[[max_index]][[5]] ## pour pi2
ggplot(as.data.frame(resultats[[max_index]][[5]]), aes(x=seq_along(resultats[[max_index]][[5]]), y=resultats[[max_index]][[5]])) +
  geom_point(color = "#21ADE5", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab(expression(pi[2]))

## pi3
resultats[[max_index]][[6]] ## pour pi3
ggplot(as.data.frame(resultats[[max_index]][[6]]), aes(x=seq_along(resultats[[max_index]][[6]]), y=resultats[[max_index]][[6]])) +
  geom_point(color = "#21ADE5", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab(expression(pi[3]))


## plot des proportions ensemble
super <- data.frame(pi1 = resultats[[max_index]][[4]], pi2 = resultats[[max_index]][[5]], pi3 = resultats[[max_index]][[6]])
super <- melt(super)
iteration <- c(1:20)
super <- cbind(super, iteration)
propplot <- ggplot(super, aes(x = iteration, y = value, colour = variable)) +
  geom_point(aes(x=iteration, y=value), size = 2.5) +
  #geom_point(aes(x=1, y=pi3[1]), color = "black", size = 1.5) +
  theme_bw() +
  xlab("Itération") +
  ylab("Proportion des groupes") +
  scale_colour_manual(values=c("hotpink1", "red", "darkgreen"), labels = c("Groupe 1", "Groupe 2", "Groupe 3")) +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="#909090")) +
  theme(legend.title = element_blank())

ggsave(plot=propplot, filename="propplot.png", width=8, height=6)



################################################################
##------------------------THETA HAT---------------------------##
################################################################

## On récupère le theta_hat final pour pouvoir classer les individus
mu_hat <- tail(resultats[[max_index]][[2]], 1)
sigma_hat <- tail(resultats[[max_index]][[3]], 1)
pi1_hat <- tail(resultats[[max_index]][[4]], 1)
pi2_hat <- tail(resultats[[max_index]][[5]], 1)
pi3_hat <- tail(resultats[[max_index]][[6]], 1)
alpha_hat <- tail(resultats[[max_index]][[7]], 1)

################################################################
##------------------------GRAPHIQUE---------------------------##
################################################################

p <- function(x){
  (pi2_hat/(pi2_hat+pi3_hat)*dcauchy(x, location = mu_hat, scale = sigma_hat) + pi3_hat/(pi2_hat+pi3_hat)*dcauchy(x, location = 2*mu_hat, scale = sqrt(2)*sigma_hat))
}

p2 <- function(x){
  pi2_hat/(pi2_hat+pi3_hat)*dcauchy(x, location = mu_hat, scale = sigma_hat)
}

p3 <- function(x){
  pi3_hat/(pi2_hat+pi3_hat)*dcauchy(x, location = 2*mu_hat, scale = sqrt(2)*sigma_hat)
}

ggplot(sum_gene_df %>% filter(nb_transcrit > 300)%>%filter(nb_transcrit < 20000), aes(x = nb_transcrit))+
  geom_histogram(aes(y = ..density.., color = "Données"), position = "identity", bins = 100, fill = "#75D7FF")+
  stat_function(fun = p2, aes(color = 'Deuxième groupe'), size = 1)+
  stat_function(fun = p3, aes(color = 'Troisième groupe'), size = 1)+
  stat_function(fun = p, aes(color = 'Mélange'), size = 1.1)+
  scale_color_manual(name = "Distributions",
                     values = c("Deuxième groupe" = "red",
                                "Troisième groupe" = "darkgreen",
                                "Mélange" = 'black',
                                "Données" =  "#65C6ED"))+
  theme(legend.background = element_rect(fill="white",
                                         size=1, linetype="solid",
                                         colour ="#909090"),
        legend.position = c(0.8, 0.9),
        legend.direction = "vertical",
        legend.title = element_blank())+
  xlab("Nombre de transcrits") +
  ylab("Densité")+
  theme_bw()

p_tout <- function(x){
  (pi1_hat*dpareto(x, k = alpha_hat, xmin = 1) + pi2_hat*dcauchy(x, location = mu_hat, scale = sigma_hat) + pi3_hat*dcauchy(x, location = 2*mu_hat, scale = sqrt(2)*sigma_hat))
}

ggplot(sum_gene_df %>% filter(nb_transcrit < 50), aes(x = nb_transcrit))+
  geom_histogram(aes(y = ..density.., color = "Données"), position = "identity", bins = 250, fill = "#75D7FF")+
  stat_function(fun = p_tout, aes(color = 'Mélange'), size = 0.7)+
  scale_color_manual(name = "Distributions",
                     values = c("Deuxième groupe" = "red",
                                "Troisième groupe" = "darkgreen",
                                "Mélange" = 'black',
                                "Données" =  "#65C6ED"))+
  theme(legend.background = element_rect(fill="white",
                                         size=1, linetype="solid",
                                         colour ="#909090"),
        legend.position = c(0.8, 0.9),
        legend.direction = "vertical",
        legend.title = element_blank())+
  xlab("Nombre de transcrits") +
  ylab("Densité")+
  theme_bw()

       
pop_classe_avant <- count(sum_gene_df, nb_transcrit)
pop_classe_avant$proportion <- (pop_classe_avant$n / sum(pop_classe_avant$n))
proportion_50 <- pop_classe_avant %>% filter(nb_transcrit < 25)

plot_pareto <- ggplot(proportion_50, aes(x = nb_transcrit)) +
  geom_col(aes(x = nb_transcrit, y = proportion),just = 0.5, width = 0.99, color = "#65C6ED", fill = "#75D7FF") +
  theme_bw()+
  stat_function(aes(color= "Densité de la loi de Pareto"),size = 1, fun = function(x) {(dpareto(x, k = alpha_hat, xmin = 1)*pi1_hat)})+
  ylab("Proportion") +
  ylim(0,0.9) +
  xlab("Nombre de transcrits") +
  scale_color_manual(values = c("Densité de la loi de Pareto" = "hotpink"))+
  theme(legend.position = "top",
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="#909090")) +
  theme(legend.position = c(0.8, 0.94),
        legend.direction = "horizontal")+
  theme(legend.title = element_blank())
  #labs(title = "Fonction de masse de la loi de Pareto superposées aux données réelles") +
  #theme(plot.title = element_text(face = "bold",size = 13, hjust = 0, vjust = 0))

plot_pareto

cauchys_pareto <- ggplot(pop_classe_avant %>% filter(nb_transcrit > 50)%>%filter(nb_transcrit < 30000), aes(x = nb_transcrit)) +
  geom_col(aes(x = nb_transcrit, y = proportion), just = 0.5, width = 0.1, color = "#65C6ED", fill = "#75D7FF") +
  theme_bw()+
  stat_function(aes(color= "distribution 1"),size = 1, fun = function(x) {(dcauchy(x, location = mu_hat, scale = sigma_hat)*pi2_hat)*6})+
  stat_function(aes(color= "distribution 2"), size = 1, fun = function(x) {(dcauchy(x, location = 2*mu_hat, scale = sqrt(2)*sigma_hat)*pi3_hat)*6})+
  scale_color_manual(name = "Distributions",
                     values = c("distribution 1" = "red",
                                "distribution 2" = "darkgreen")) +
  ylab("Proportion") +
  xlab("Nombre de transcrits")

cauchys_pareto


ggplot(sum_gene_df %>% filter(nb_transcrit > 500)%>%filter(nb_transcrit < 50000), aes(x = nb_transcrit)) +
  geom_histogram(aes(x = nb_transcrit), binwidth = 100, color = "#65C6ED", fill = "#75D7FF") +
  theme_bw()+
  stat_function(aes(color= "Première normale"),size = 1, fun = function(x) {(dcauchy(x, location = mu_hat, scale = sigma_hat)*pi2_hat*34500000)})+
  stat_function(aes(color= "Seconde normale"), size = 1, fun = function(x) {(dcauchy(x, location = 2*mu_hat, scale = sqrt(2)*sigma_hat)*pi3_hat*34500000)})+
  scale_color_manual(name = "Distributions",
                     values = c("Première normale" = "red",
                                "Seconde normale" = "darkgreen")) +
  ylab("Effectif") +
  theme(legend.position = "top",
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="#909090"))+
  xlab("Nombre de transcrits") +
  theme(legend.position = c(0.8, 0.9),
        legend.direction = "vertical")+
  theme(legend.title = element_blank())



################################################################
##----------------------Classification------------------------##
################################################################

## Les probabilités sa caculent par le k qui maximise le tik(theta)
## Il faut donc calculer les tik(theta) pour les individus et tous les groupes
logphi1 <- dpareto(data, k = alpha_hat, xmin = 1, log = TRUE)
logphi2 <- logspace.sub(pcauchy(data+1/2, location = mu_hat, scale = sigma_hat, log.p = TRUE), pcauchy(data-1/2, location = mu_hat, scale = sigma_hat, log.p = TRUE))
logphi3 <- logspace.sub(pcauchy(data+1/2, location = 2*mu_hat, scale = sqrt(2)*sigma_hat, log.p = TRUE), pcauchy(data-1/2, location = 2*mu_hat, scale = sqrt(2)*sigma_hat, log.p = TRUE))


ln1 <- log(pi1_hat) + logphi1
ln2 <- log(pi2_hat) + logphi2
ln3 <- log(pi3_hat) + logphi3

somme_phi_pondere <- pi1_hat*exp(logphi1) + pi2_hat*exp(logphi2) + pi3_hat*exp(logphi3)

t1 <- pi1_hat*exp(logphi1) / somme_phi_pondere
t2 <- pi2_hat*exp(logphi2) / somme_phi_pondere
t3 <- pi3_hat*exp(logphi3) / somme_phi_pondere


## On regroupe les probas de chaque indivdus pour chaque groupe
proba <- cbind(t1, t2, t3)

## On créé la liste des groupes attribués à chaque individus
groupe <- case_when(proba[,2]>0.99 ~ 2,
                    proba[,1]>0.8 ~ 1,
                    proba[,3]>0.2 ~ 3,
                    TRUE ~ 1)

## On se donne un dataframe constitué des individus et de leur groupe associé
data_classee <- as.data.frame(cbind(data, groupe))

## On a un summary par groupe
tapply(data_classee$data, data_classee$groupe, sd)

which(data_classee$groupe==2)

## On compte le nombre d'individu par groupe, et leur proportion (on retrouve les proportions du modeles)
population_par_classe <- count(data_classee, groupe)
population_par_classe$proportion <- round((population_par_classe$n / length(data))*100, 3)
population_par_classe

## On peut compter le nombre d'individu par groupe et leur proportion avant le modèle
pop_classe_avant <- count(sum_gene_df, classe)
pop_classe_avant$proportion <- round((pop_classe_avant$n / length(data))*100, 3)
pop_classe_avant


## On compte le nombre d'occurence de chaque individus et leur proportion par rapport à la population totale
## Piste pour amélioration, est ce que on peut retirer les premières gouttelettes ?
population_par_sumgene <- count(data_classee, data)
population_par_sumgene$proportion <- round((population_par_sumgene$n / length(data))*100, 3)
population_par_sumgene

plot_classe <- ggplot(data = data_classee) +
  geom_boxplot(aes(x = factor(groupe), y = data)) +
  scale_x_discrete(breaks = 1:3, labels = c("Groupe 1", "Groupe 2", "Groupe 3")) +
  #scale_y_continuous (trans='log10') +
  labs(title = "Etude de la performance des algorithmes d'optimisation",
       x = "Taille de l'échantillon",
       y = "Temps de convergence (échelle logarithmique)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  theme_bw()
plot_classe



################################################################
##------------------------Recherche---------------------------##
################################################################



## On peut générer un échantillon de notre loi de mélange et tracer son histogramme
taille_echantillon <- 350000
echantillon_hat <- c(rpois(taille_echantillon*pi1_hat, lambda_hat+3), round(rnorm(taille_echantillon*pi2_hat, mu_hat, sigma_hat)), round(rnorm(taille_echantillon*pi3_hat, 2*mu_hat, sqrt(2)*sigma_hat)))
echantillon_hat <- as.data.frame(cbind(echantillon_hat, c(rep(1, taille_echantillon*pi1_hat), rep(2, taille_echantillon*pi2_hat), rep(3, taille_echantillon*pi3_hat))))
colnames(echantillon_hat) <- c("echantillon", "groupe")
echantillon_hat <- echantillon_hat[which(echantillon_hat$echantillon > 0),]
echantillon_hat <- echantillon_hat[sample(nrow(echantillon_hat)),]

ggplot(echantillon_hat, aes(x = seq_along(echantillon), y = echantillon, color = as.factor(groupe))) +
  geom_point(shape = '.', alpha = 0.4) +
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

## Groupe 1
echantillon_hat_1 <- echantillon_hat[which(echantillon_hat$groupe == 1),]
echantillon_hat_1 <- subset(echantillon_hat_1, select = -groupe)

data_classee_1 <- data_classee[which(data_classee$groupe == 1),]
data_classee_1 <- subset(data_classee_1, select = -groupe)

DF1 <- data.frame(variable=rep(c('Echantillon issu de la loi de mélange', 'Données réelles'), each=300000), value=c(echantillon_hat_1$echantillon[1:300000] , data_classee_1$data[1:300000]))

ggplot(DF1) + 
  stat_ecdf(aes(value, color=variable), size = 1) +
  scale_color_manual(values=c("red", "#66CCFF")) +
  theme_bw()

tapply(DF1$value, DF1$variable, summary)


## Groupe 2
echantillon_hat_2 <- echantillon_hat[which(echantillon_hat$groupe == 2),]
echantillon_hat_2 <- subset(echantillon_hat_2, select = -groupe)

data_classee_2 <- data_classee[which(data_classee$groupe == 2),]
data_classee_2 <- subset(data_classee_2, select = -groupe)

DF2 <- data.frame(variable=rep(c('Echantillon issu de la loi de mélange', 'Données réelles'), each=300000), value=c(echantillon_hat_2$echantillon[1:300000] , data_classee_2$data[1:300000]))

ggplot(DF2) + 
  stat_ecdf(aes(value, color=variable), size = 1) +
  scale_color_manual(values=c("red", "#66CCFF")) +
  theme_bw()

tapply(DF2$value, DF2$variable, summary)

## Groupe 3
echantillon_hat_3 <- echantillon_hat[which(echantillon_hat$groupe == 3),]
echantillon_hat_3 <- subset(echantillon_hat_3, select = -groupe)

data_classee_3 <- data_classee[which(data_classee$groupe == 3),]
data_classee_3 <- subset(data_classee_3, select = -groupe)

DF3 <- data.frame(variable=rep(c('Echantillon issu de la loi de mélange', 'Données réelles'), each=300000), value=c(echantillon_hat_3$echantillon[1:300000] , data_classee_3$data[1:300000]))

ggplot(DF3) + 
  stat_ecdf(aes(value, color=variable), size = 1) +
  scale_color_manual(values=c("red", "#66CCFF")) +
  theme_bw()

tapply(DF3$value, DF3$variable, summary)

###

DF <- data.frame(variable=rep(c('Echantillon issu de la loi de mélange', 'Données réelles'), each=length(echantillon_hat)), value=c(echantillon_hat$echantillon , data[1:length(echantillon_hat)]))

ggplot(DF) + 
  stat_ecdf(aes(value, color=variable), size = 1) +
  scale_color_manual(values=c("red", "#66CCFF")) +
  theme_bw()



ggplot(as.data.frame(echantillon_hat), aes(x = echantillon_hat)) + 
  stat_ecdf(geom = "point",  color = "red") +
  labs(title="Cumulative Density Function \nfor Poisson",
       y = "F(x)", x="x") +
  theme_bw()

ggplot(as.data.frame(data), aes(x = data)) + 
  stat_ecdf(geom = "point",  color = "red") +
  labs(title="Cumulative Density Function \nfor Poisson",
       y = "F(x)", x="x") +
  theme_bw()




