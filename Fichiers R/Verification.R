library(DPQ) #logspace.sub
library(nloptr) #optimisation numérique
library(doParallel) #créer les clusters
library(foreach) #permet le parallélisme
library(reshape2) #formate les données pour ggplot
library(ggplot2) #créer les boxplots
library(tictoc) #compte le temps

################################################################################
##------------------------------INSTRUCTIONS----------------------------------##
################################################################################
# Ce fichier contient la vérification que l'algorithme EM converge bien vers les
# bon paramètres, pour un modèle de mélange contenant une loi de Poisson 
# et deux lois normales.

# Dans ce fichier les résultats n'ont pas été sauvegarder. 

# Les paramètres qui influencent le temps d'execution sont la taille 
# d'échantillon et le nombre de random start.
# Le temps d'execution est d'autant plus rapide que le nombre de coeur est 
# proche du nombre de répétiton. Mais il n'y a aucun intérêt à avoir un nombre 
# de coeur supérieur au nombre de répétition.

# L'analyse des résultats se fait ensuite à la main par manque de temps.



################################################################################
##-------------------------------FONCTIONS------------------------------------##
################################################################################
rm(list=ls())
set.seed(1664)

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
##-------------------------------ALGORITHME-----------------------------------##
################################################################################

## Initialisation du parallélisme
nbre_coeurs_disponibles = detectCores()
nbre_coeurs_disponibles
nbre_de_coeurs_voulu <- 20  ## Il faut que ca soit inférieur au nombre de coeurs disponibles
nbre_coeurs_voulu <- makeCluster(nbre_de_coeurs_voulu)
registerDoParallel(nbre_coeurs_voulu)

## Initialisation des paramètres des lois cible
lambda_cible <- 50
mu_cible <- 100
sigma_cible <- 10

taille_echantillon <- c(10^4)
nbre_repetition <- 20
algorithmes <- c("NLOPT_LN_BOBYQA")  
nbre_random_start <- 10 #### 2x le nombre de paramètres.
nbre_parametre_interet <- 9 #### log-vraisemblance complétée, lambda, mu, sigma, pi1, pi2, pi3

#!!# On ne s'intéresse ici que aux resultats finaux. On ne s'intéresse pas à l'évolution de la log-vraisemblance complétée ni à l'évolution de theta.
#!!# En revanche, ces aspects sont fondamentals dans le rapport et dans l'exécution de l'algorithme sur nos données réelles.

resultats_finaux <- array(data = NA, c(nbre_parametre_interet, length(taille_echantillon), nbre_repetition, length(algorithmes), nbre_random_start))
# resultats_finaux[1:9,i,j,k,l]

### Des commentaires sur le parallélisme ###
### La boucle foreach doit se trouver en première. Ca permet de réserver les coeurs.
### Il faut inclure les packages dans la boucle car elle créer des sessions R différentes. Et pour faire tourner notre code on a besoin de ces packages. C'est aussi pour ça qu'il faut changer le setseed et initialiser le resultats_temp dans la boucle pour que chaque worker ait ses infos.

# On commence à compter le temps
tic("Vérification")

all_results_temp <-foreach (j=1:nbre_repetition, .packages = c("DPQ", "nloptr")) %dopar% {
  set.seed(j)
  resultats_temp <- array(data = NA, c(nbre_parametre_interet, length(taille_echantillon), length(algorithmes), nbre_random_start))
  
  for (i in 1:length(taille_echantillon)) {
    taille_echantillon_en_cours <- taille_echantillon[i]
    data <- c(rpois(taille_echantillon_en_cours*0.6, lambda_cible), round(rnorm(taille_echantillon_en_cours*0.1, mu_cible, sigma_cible)), round(rnorm(taille_echantillon_en_cours*0.3, 2*mu_cible, sqrt(2)*sigma_cible)))
    
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

# On arrête de compter le temps
toc()

# On stop l'utilisation de plusieurs cluster.
stopCluster(nbre_coeurs_voulu)

################################################################################
##--------------------------------ANALYSE-------------------------------------##
################################################################################
## On récupère un tableau juste avec la log-vraisemblance complétée
loglikelihood <- resultats_finaux[1, , , , ]
loglikelihood

lambda_all <- resultats_finaux[4, , , , ]
lambda_all
lambda_liste <- c()
for (i in 1:length(lambda_all)) {
  if (!is.na(lambda_all[i])) { if (lambda_all[i] < 51 & lambda_all[i] > 48) {lambda_liste <- c(lambda_liste, lambda_all[i])} }
}
mean(lambda_liste)
sd(lambda_liste)

mu_all <- resultats_finaux[5, , , , ]
mu_all
mu_liste <- c()
for (i in 1:length(mu_all)) {
  if (!is.na(mu_all[i])) { if (mu_all[i] < 101 & mu_all[i] > 99) {mu_liste <- c(mu_liste, mu_all[i])} }
}
mean(mu_liste)
sd(mu_liste)

sigma_all <- resultats_finaux[6, , , , ]
sigma_all
sigma_liste <- c()
for (i in 1:length(sigma_all)) {
  if (!is.na(sigma_all[i])) { if (sigma_all[i] < 11 & sigma_all[i] > 9.5) {sigma_liste <- c(sigma_liste, sigma_all[i])} }
}
mean(sigma_liste)
sd(sigma_liste)

pi1_all <- resultats_finaux[7, , , , ]
pi1_all
pi1_liste <- c()
for (i in 1:length(pi1_all)) {
  if (!is.na(pi1_all[i])) { if (pi1_all[i] < 0.62 & pi1_all[i] > 0.58) {pi1_liste <- c(pi1_liste, pi1_all[i])} }
}
mean(pi1_liste)
sd(pi1_liste)

pi2_all <- resultats_finaux[8, , , , ]
pi2_all
pi2_liste <- c()
for (i in 1:length(pi2_all)) {
  if (!is.na(pi2_all[i])) { if (pi2_all[i] < 0.11 & pi2_all[i] > 0.09) {pi2_liste <- c(pi2_liste, pi2_all[i])} }
}
mean(pi2_liste)
sd(pi2_liste)

pi3_all <- resultats_finaux[9, , , , ]
pi3_all
pi3_liste <- c()
for (i in 1:length(pi3_all)) {
  if (!is.na(pi3_all[i])) { if (pi3_all[i] < 0.32 & pi3_all[i] > 0.28) {pi3_liste <- c(pi3_liste, pi3_all[i])} }
}
mean(pi3_liste)
sd(pi3_liste)
