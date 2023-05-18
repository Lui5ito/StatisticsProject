library(DPQ) #logspace.sub
library(nloptr) #nolptr
library(doParallel) #créer les clusters
library(foreach) #permet le parallélisme
library(reshape2) #formate les données pour ggplot
library(ggplot2) #créer les boxplots

################################################################################
##------------------------------INSTRUCTIONS----------------------------------##
################################################################################
# Ce fichier contient la comparaison entre différents algortithmes
# d'optimisation, pour un modèle de mélange avec une loi de Poisson tronquée en
# zéro et deux lois normales.

# La partie 'ALGORITHME' mets du temps à tourner. Pour explorer les résultats,
# il suffit de lancées la parties 'INITIALISATION', puis la partie 'RESULTATS'.

rm(list=ls())
set.seed(1664)
################################################################################
##-------------------------------FONCTIONS------------------------------------##
################################################################################


logvraissemblance <- function(params){
  logphi1 <- dpois(data, lambda_r, log = TRUE)
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
##-----------------------------INITIALISATION---------------------------------##
################################################################################
## Initialisation des paramètres des lois cible
lambda_cible <- 3
mu_cible <- 5050
sigma_cible <- 1010


taille_echantillon <- c(10^2, 10^3, 10^4, 10^5)
nbre_repetition <- 20
algorithmes <- c("NLOPT_LN_NELDERMEAD", "NLOPT_LN_COBYLA", "NLOPT_LN_BOBYQA")
nbre_random_start <- 10
nbre_parametre_interet <- 9 #### log-vraisemblance complétée, temps, iterations, lambda, mu, sigma, pi1, pi2, pi3

################################################################################
##-------------------------------ALGORITHME-----------------------------------##
################################################################################

## Initialisation du parallélisme
nbre_coeurs_disponibles = detectCores()
nbre_coeurs_disponibles
nbre_de_coeurs_voulu <- 20  ## Il faut que ca soit inférieur au nombre de coeurs disponibles
nbre_coeurs_voulu <- makeCluster(nbre_de_coeurs_voulu)
registerDoParallel(nbre_coeurs_voulu)


#!!# On ne s'intéresse ici que aux resultats finaux. On ne s'intéresse pas à l'évolution de la log-vraisemblance complétée ni à l'évolution de theta.
#!!# En revanche, ces aspects sont fondamentaux dans le rapport et dans l'exécution de l'algorithme sur nos données réelles.

resultats_finaux <- array(data = NA, c(nbre_parametre_interet, length(taille_echantillon), nbre_repetition, length(algorithmes), nbre_random_start))
# resultats_finaux[1:9,i,j,k,l]

### Des commentaires sur le parallélisme ###
### La boucle foreach doit se trouver en première. Ca permet de réserver les coeurs.
### Il faut inclure les packages dans la boucle car elle créer des sessions R différentes. Et pour faire tourner notre code on a besoin de ces packages. C'est aussi pour ça qu'il faut changer le setseed et initialiser le resultats_temp dans la boucle pour que chaque worker ait ses infos.

all_results_temp <-foreach (j=1:nbre_repetition, .packages = c("DPQ", "nloptr")) %dopar% {
  set.seed(j)
  resultats_temp <- array(data = NA, c(nbre_parametre_interet, length(taille_echantillon), length(algorithmes), nbre_random_start))
  
  for (i in 1:length(taille_echantillon)) {
    taille_echantillon_en_cours <- taille_echantillon[i]
    data <- c(rpois(taille_echantillon_en_cours/2, lambda_cible), round(rnorm(taille_echantillon_en_cours/4, mu_cible, sigma_cible)), round(rnorm(taille_echantillon_en_cours/4, 2*mu_cible, sqrt(2)*sigma_cible)))
    
    for (k in 1:length(algorithmes)){
      algo_en_cours <- algorithmes[k]
      print(algo_en_cours)
      
      ## Cette boucle random start EST l'algorithme EM
      for (l in 1:nbre_random_start){
        
        n <- length(data) ## Ici n vaut i mais j'ai déjà écrit n partout dans la suite...
        
        ## Initialisation aléatoire des paramètres
        lambda_r <- data[sample(1:n, 1)]
        mu_r <- data[sample(1:n, 1)]
        sigma_r <- 100*sd(data)
        pi1_r <- 1/3
        pi2_r <- 1/3
        pi3_r <- 1/3
        Lvc <- c()
        temps_nlopt <- c()
        iteration_nlopt <- c()
        
        t <- try(repeat{
          
          ## On a bseoin de ces calculs pour lambda_r
          logphi1 <- dpois(data, lambda_r, log = TRUE)
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
          
          ## On peut maintenant mettre à jour lambda_r
          lambda_r <- sum(t1*data)/T1
          
          ## On doit calculer la log-vraisemblance à chaque itération puisque c'est notre porte de sortie de la boucle repeat
          LV_r <- lv1+lv2+lv3
          Lvc <- c(Lvc, LV_r)
          
          ## On optimise avec l'algorithme k
          temps_ce_tour <- system.time(res_nlopt <- nloptr(x0 = c(mu_r, sigma_r), eval_f = logvraissemblance_pour_nloptr, opts = list(algorithm = algo_en_cours, maxeval = 10000, tol_rel=1e-15, xtol_abs=1e-15), lb = c(1, 1)))
          temps_nlopt <- c(temps_nlopt, temps_ce_tour[3])
          iteration_nlopt <- c(iteration_nlopt, res_nlopt$iteration)
          
          ## Mises à jours des paramètres
          mu_r <- res_nlopt$solution[1]
          sigma_r <- res_nlopt$solution[2]
          pi1_r <- T1/n
          pi2_r <- T2/n
          pi3_r <- T3/n
          
          if (length(Lvc)>2 & abs((tail(Lvc, 1) - tail(Lvc, 2)[1])) < 0.1) { break }
        })
        
        
        ## On vérifie qu'il n'y a pas eu une seule erreur de nlopt dans l'algorithme EM
        if(!(inherits(t, "try-error"))){
          ## On sauvegarde alors tous les paramètres d'intérêts finaux de notre algorithme EM
          resultats_temp[1, i, k, l] <- logvraissemblance(c(mu_r, sigma_r)) ## log-vraisemblance complétée
          resultats_temp[2, i, k, l] <- sum(temps_nlopt) ## On prend la moyenne de tous les nlopt utilisés lors de cet algorithme EM
          resultats_temp[3, i, k, l] <- sum(iteration_nlopt) ## On doit encore prendre la moyenne du nombre d'itération de nlopt lors de l'algorithme EM
          resultats_temp[4, i, k, l] <- lambda_r
          resultats_temp[5, i, k, l] <- mu_r
          resultats_temp[6, i, k, l] <- sigma_r
          resultats_temp[7, i, k, l] <- pi1_r
          resultats_temp[8, i, k, l] <- pi2_r
          resultats_temp[9, i, k, l] <- pi3_r
        }
        
      }
    }
  }
  return(resultats_temp)
}


## On remet tous les calculs parallèle ensemble.
for (j in 1:nbre_repetition){
  resultats_finaux[,,j,,] <- all_results_temp[[j]]
}

## On stop l'utilisation de plusieurs cluster.
stopCluster(nbre_coeurs_voulu)

## On sauvegarde les résultats.
saveRDS(resultats_finaux, file = "resultats_comparaison.RData")

################################################################################
##--------------------------------RESULTATS-----------------------------------##
################################################################################

resultats_finaux <- readRDS("resultats_comparaison.RData")

summary(resultats_finaux)

## On récupère un tableau juste avec la log-vraisemblance complétée
loglikelihood <- resultats_finaux[1, , , , ]

max_loglikelihood <- apply(X = loglikelihood, MARGIN = c(1,2,3), max, na.rm = TRUE)
max_max_loglikelihood <- apply(X = max_loglikelihood, MARGIN = c(1,3), max, na.rm = TRUE)

min_loglikelihood <- apply(X = loglikelihood, MARGIN = c(1,2,3), min, na.rm = TRUE)
min_min_loglikelihood <- apply(X = min_loglikelihood, MARGIN = c(1,3), min, na.rm = TRUE)

## On récupère la log-vraisemblance compltée moyenne sur les random start
mean_loglikelihood <- apply(X = loglikelihood, MARGIN = c(1,2,3), mean, na.rm = TRUE)
mean_mean_loglikelihood <- apply(X = mean_loglikelihood, MARGIN = c(1,3), mean, na.rm = TRUE)

## On récupère aussi les écart-type
sd_loglikelihood <- apply(X = loglikelihood, MARGIN = c(1,2,3), sd, na.rm = TRUE)
sd_sd_loglikelihood <- apply(X = sd_loglikelihood, MARGIN = c(1,3), sd, na.rm = TRUE)


## On récupère un tableau avec juste le temps
time <- resultats_finaux[2, , , , ]

## Table formatée pour le ggplot
for_boxplot <- melt(time, varnames = c("taille_echantillon", "nbre_repetition", "algorithmes", "nbre_random_start"))

## Le ggplot
plot_temps <- ggplot(data = for_boxplot) +
  geom_boxplot(aes(x = factor(taille_echantillon), y = value, fill = factor(algorithmes))) +
  scale_x_discrete(breaks = 1:length(taille_echantillon), labels = taille_echantillon) +
  scale_y_continuous (trans='log10') +
  labs(title = "Etude de la performance des algorithmes d'optimisation",
       x = "Taille de l'échantillon",
       y = "Temps de convergence (échelle logarithmique)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  labs(fill = "") +
  theme(legend.position = "right") +
  scale_fill_manual(
    values = c("#93dbff", "#66CCFF", "#51a3cc"),
    name = "Algorithmes", 
    labels = c("NelderMead", "Cobyla", "Bobyqa")
  ) +
  theme_bw()

plot_temps

