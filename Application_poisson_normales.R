library(dplyr) #nécessaire pour l'importation des données
library(DPQ) #logspace.sub
library(nloptr) #optimisation numérique
library(doParallel) #créer les clusters
library(foreach) #permet le parallélisme
library(tictoc) #compte le temps
library(ggplot2) #graphs

################################################################################
##------------------------------INSTRUCTIONS----------------------------------##
################################################################################
# Ce fichier contient l'algorithme EM pour un modèle de mélange contenant 
# une loi de Poisson tronquée en zéro et deux lois normales.

# La partie 'ALGORITHME' mets du temps à tourner. Pour explorer les résultats,
# il suffit de lancées la parties 'DONNEES', puis la partie 'RESULTATS'.

################################################################################
##--------------------------------DONNEES-------------------------------------##
################################################################################
rm(list=ls())
set.seed(1664)

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
  logphi1 <- dpois(data, params[3], log = TRUE) - logspace.sub(0, -params[3])
  logphi2 <- logspace.sub(pnorm(data+1/2, mean = params[1], sd = params[2], log.p = TRUE), pnorm(data-1/2, mean = params[1], sd = params[2], log.p = TRUE))
  logphi3 <- logspace.sub(pnorm(data+1/2, mean = 2*params[1], sd = sqrt(2)*params[2], log.p = TRUE), pnorm(data-1/2, mean = 2*params[1], sd = sqrt(2)*params[2], log.p = TRUE))
  
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

nbre_random_start <- 100
nbre_parametre_interet <- 7 ## log-vraisemblance complétée, lambda, mu, sigma, pi1, pi2, pi3


## Initialisation de tous ce qui est statique
n <- length(data)
n_sup_1000 <- length(data_sup_1000)

## On commence à compter le temps
tic("Application")

## Cette boucle random start EST l'algorithme EM
## resultats est une liste contenant chaque random start. 
## Un random start est une liste contenant des vecteurs qui contiennent eux même les itérations de chaque paramètre d'interêt.

resultats <- foreach (i=1:nbre_random_start, .packages=c("DPQ", "nloptr")) %dopar% {
  ## On setseed à chaque tour pour changer les initilisations.
  set.seed(i)
  
  ## On initialise la liste qui va contenir le résultat du random start de l'algorithme 
  ma_liste <- NULL
  
  ## (ré-)Initialisation des suivis des paramètres d'interets
  suite_lambda <- c()
  suite_mu <- c()
  suite_sigma <- c()
  suite_pi1 <- c()
  suite_pi3 <- c()
  suite_pi2 <- c()
  Lvc <- c()
  
  ## Initialisation aléatoire des paramètres
  lambda_r <- data[sample(1:n, 1)]
  mu_r <- data_sup_1000[sample(1:n_sup_1000, 1)] 
  sigma_r <- data_sup_1000[sample(1:n_sup_1000, 1)]
  pi1_r <- 1/3
  pi2_r <- 1/3
  pi3_r <- 1/3
  
  ## Première instance des paramètres d'interets
  suite_lambda <- c(lambda_r)
  suite_sigma <- c(sigma_r)
  suite_mu <- c(mu_r)
  suite_pi1 <- c(pi1_r)
  suite_pi2 <- c(pi2_r)
  suite_pi3 <- c(pi3_r)
  
  t <- try(repeat{
    
    logphi1 <- dpois(data, lambda_r, log = TRUE) - logspace.sub(0, -lambda_r)
    logphi2 <- logspace.sub(pnorm(data+1/2, mean = mu_r, sd = sigma_r, log.p = TRUE), pnorm(data-1/2, mean = mu_r, sd = sigma_r, log.p = TRUE))
    logphi3 <- logspace.sub(pnorm(data+1/2, mean = 2*mu_r, sd = sqrt(2)*sigma_r, log.p = TRUE), pnorm(data-1/2, mean = 2*mu_r, sd = sqrt(2)*sigma_r, log.p = TRUE))
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
    res_nlopt <- nloptr(x0 = c(mu_r, sigma_r, lambda_r), eval_f = logvraissemblance_pour_nloptr, opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 10000, tol_rel=1e-15, xtol_abs=1e-15), lb = c(1, 1, 0.1))
    
    ## Mises à jours des paramètres
    mu_r <- res_nlopt$solution[1]
    sigma_r <- res_nlopt$solution[2]
    lambda_r <- res_nlopt$solution[3]
    pi1_r <- T1/n
    pi2_r <- T2/n
    pi3_r <- T3/n
    
    ## Sauvegarde des historiques des paramètres d'intérêts
    suite_lambda <- c(suite_lambda, lambda_r)
    suite_mu <- c(suite_mu, mu_r)
    suite_sigma <- c(suite_sigma, sigma_r)
    suite_pi1 <- c(suite_pi1, pi1_r)
    suite_pi2 <- c(suite_pi2, pi2_r)
    suite_pi3 <- c(suite_pi3, pi3_r)
    
    ## Condition d'arrêt : la log-vraisemblance complétée n'a pas augmentée de plus de 0.1
    if (length(Lvc)>2 & abs((tail(Lvc, 1) - tail(Lvc, 2)[1])) < 0.1) { break }
  })
  
  
  ## On vérifie qu'il n'y a pas eu une seule erreur de nlopt dans le repeat
  if(!(inherits(t, "try-error"))){
    
    ## On sauvegarde alors tous les paramètres d'intérêts finaux de notre algorithme EM
    Lvc <- c(Lvc, logvraissemblance(c(mu_r, sigma_r, lambda_r))) #Pour avoir autant de points que les paramètres
    ma_liste <- list(logvraisemblance  = Lvc, lambda = suite_lambda, mu = suite_mu, sigma = suite_sigma, pi1 = suite_pi1, pi2 = suite_pi2, pi3 = suite_pi3)
  }
  return (ma_liste)
}

## On arrête de compter le temps
toc()

## On arrête le parallélisme
stopCluster(nbre_coeurs_voulu)

## On sauvagarde le résultat de l'algorithme
saveRDS(resultats, file = "resultats_poisson_normales.RData")

################################################################################
##--------------------------------RESULTATS-----------------------------------##
################################################################################

resultats <- readRDS("resultats_poisson_normales.RData")

## On veut maintenant récupérer le meilleur des random starts, celui qui a la log-vraisemblance complétée la plus élevée.
max_index <- 1
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

## lambda
resultats[[max_index]][[2]]
iteration_lambda <- ggplot(as.data.frame(resultats[[max_index]][[2]]), aes(x=seq_along(resultats[[max_index]][[2]]), y=resultats[[max_index]][[2]])) +
  geom_point(color = "#21ADE5", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab(expression(lambda))
ggsave(plot=iteration_lambda, filename="iteration_lambda.png", width=8, height=5)

## mu
resultats[[max_index]][[3]]
iteration_mu <- ggplot(as.data.frame(resultats[[max_index]][[3]]), aes(x=seq_along(resultats[[max_index]][[3]]), y=resultats[[max_index]][[3]])) +
  geom_point(color = "#21ADE5", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab(expression(mu))
ggsave(plot=iteration_mu, filename="iteration_mu.png", width=8, height=5)

## sigma
resultats[[max_index]][[4]]
iteration_sigma <- ggplot(as.data.frame(resultats[[max_index]][[4]]), aes(x=seq_along(resultats[[max_index]][[4]]), y=resultats[[max_index]][[4]])) +
  geom_point(color = "#21ADE5", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab(expression(sigma))
ggsave(plot=iteration_sigma, filename="iteration_sigma.png", width=8, height=5)


## pi1
resultats[[max_index]][[5]]
ggplot(as.data.frame(resultats[[max_index]][[5]]), aes(x=seq_along(resultats[[max_index]][[5]]), y=resultats[[max_index]][[5]])) +
  geom_point(color = "#21ADE5", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab(expression(pi[1]))

## pi2
resultats[[max_index]][[6]] ## pour pi2
ggplot(as.data.frame(resultats[[max_index]][[6]]), aes(x=seq_along(resultats[[max_index]][[6]]), y=resultats[[max_index]][[6]])) +
  geom_point(color = "#21ADE5", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab(expression(pi[2]))

## pi3
resultats[[max_index]][[7]] ## pour pi3
ggplot(as.data.frame(resultats[[max_index]][[7]]), aes(x=seq_along(resultats[[max_index]][[7]]), y=resultats[[max_index]][[7]])) +
  geom_point(color = "#21ADE5", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab(expression(pi[3]))


## graphique des proportions ensemble
pi_123 <- data.frame(pi1 = resultats[[max_index]][[5]], pi2 = resultats[[max_index]][[6]], pi3 = resultats[[max_index]][[7]])
pi_123 <- melt(super)
iteration <- c(1:length(pi_123[,1]))
pi_123 <- cbind(pi_123, iteration)
prop_plot <- ggplot(pi_123, aes(x = iteration, y = value, colour = variable)) +
  geom_point(aes(x=iteration, y=value), size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab("Proportion des groupes") +
  scale_colour_manual(values=c("hotpink1", "red", "darkgreen"), labels = c("Groupe 1", "Groupe 2", "Groupe 3")) +
  theme(legend.position="top", legend.box = "horizontal") +
  theme(legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="#909090")) +
  theme(legend.title = element_blank())
ggsave(plot=prop_plot, filename="prop_plot.png", width=8, height=6)


################################################################
##------------------------THETA HAT---------------------------##
################################################################

## On récupère le theta_hat final pour pouvoir classer les individus
lambda_hat <- tail(resultats[[max_index]][[2]], 1)
mu_hat <- tail(resultats[[max_index]][[3]], 1)
sigma_hat <- tail(resultats[[max_index]][[4]], 1)
pi1_hat <- tail(resultats[[max_index]][[5]], 1)
pi2_hat <- tail(resultats[[max_index]][[6]], 1)
pi3_hat <- tail(resultats[[max_index]][[7]], 1)

################################################################
##------------------------GRAPHIQUES--------------------------##
################################################################

## Graphique pour les groupes 2 et 3
melange_groupe_2_3 <- function(x){
  (pi2_hat/(pi2_hat+pi3_hat)*dnorm(x, mean = mu_hat, sd = sigma_hat) + pi3_hat/(pi2_hat+pi3_hat)*dnorm(x, mean = 2*mu_hat, sd = sqrt(2)*sigma_hat))
}

graph_groupes_2_3 <- ggplot(sum_gene_df %>% filter(nb_transcrit > 300)%>%filter(nb_transcrit < 20000), aes(x = nb_transcrit))+
  geom_histogram(aes(y = ..density.., color = "Données"), position = "identity", bins = 100, fill = "#75D7FF")+
  stat_function(fun = melange_groupe_2_3, aes(color = 'Mélange'), size = 1.1)+
  theme_bw()+
  scale_color_manual(name = "Distributions",
                     values = c("Mélange" = 'black',
                                "Données" =  "#65C6ED"))+
  theme(legend.background = element_rect(fill="white",
                                         size=1, linetype="solid",
                                         colour ="#909090"),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title = element_blank())+
  xlab("Nombre de transcrits") +
  ylab("Densité")+
  labs(title = "Densité de la loi de mélange superposée aux données réelles") +
  theme(plot.title = element_text(face = "bold",size = 13, hjust = 0, vjust = 0))

graph_groupes_2_3

## Graphique pour le groupe 1
## Ce graph est un peu différent car il fait intervenir une fonction de masse

# On calcul les proportions à la main, ainsi que les valeurs de la fonction de masse de Poisson.
pop_classe_avant <- count(sum_gene_df, nb_transcrit)
pop_classe_avant$proportion <- (pop_classe_avant$n / sum(pop_classe_avant$n))
proportion_25 <- pop_classe_avant %>% filter(nb_transcrit < 25)
proportion_25$poisson <- exp(dpois(1:24, lambda_hat, log = TRUE)  - logspace.sub(0, -lambda_hat))


graph_groupe_1 <- ggplot(proportion_25, aes(x = nb_transcrit)) +
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
  theme(legend.position = c(0.75, 0.9),
        legend.direction = "vertical")+
  theme(legend.title = element_blank())+
  labs(title = "Fonction de masse de la loi de Poisson superposée aux données réelles") +
  theme(plot.title = element_text(face = "bold",size = 13, hjust = 0, vjust = 0))

graph_groupe_1

################################################################
##----------------------Classification------------------------##
################################################################

## Les probabilités sa caculent par le k qui maximise le tik(theta)
## Il faut donc calculer les tik(theta) pour les individus et tous les groupes
logphi1 <- dpois(data, lambda_hat, log = TRUE) - logspace.sub(0, -lambda_hat)
logphi2 <- logspace.sub(pnorm(data+1/2, mean = mu_hat, sd = sigma_hat, log.p = TRUE), pnorm(data-1/2, mean = mu_hat, sd = sigma_hat, log.p = TRUE))
logphi3 <- logspace.sub(pnorm(data+1/2, mean = 2*mu_hat, sd = sqrt(2)*sigma_hat, log.p = TRUE), pnorm(data-1/2, mean = 2*mu_hat, sd = sqrt(2)*sigma_hat, log.p = TRUE))

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
groupe <- max.col(proba)

## On se donne un dataframe constitué des individus et de leur groupe associé
data_classee <- as.data.frame(cbind(data, groupe))

## On a un summary par groupe
tapply(data_classee$data, data_classee$groupe, summary)
