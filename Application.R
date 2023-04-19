library(dplyr) #nécessaire pour l'importation des données
library(DPQ) #logspace.sub
library(nloptr) #nolptr
library(doParallel) #créer les clusters
library(foreach) #permet le parallélisme
library(tictoc)

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
pourcentage_des_coeurs_voulu <- 0.75
nbre_coeurs_voulu <- makeCluster(pourcentage_des_coeurs_voulu*nbre_coeurs_disponibles)
registerDoParallel(nbre_coeurs_voulu)

nbre_random_start <- 100 #### On prend plus que 2x le nombre de paramètres car on a un échantillon dispoportionné
nbre_parametre_interet <- 7 #### log-vraisemblance complétée, lambda, mu, sigma, pi1, pi2, pi3


## Initialisation de tous ce qui est statique
n <- length(data)
n_sup_1000 <- length(data_sup_1000)

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
    ## Sauvegarde de la valeur de la log-vraisemblance complétée à cette instant
    Lvc <- c(Lvc, LV_r)
    
    ## On optimise avec l'algorithme Bobyqat
    res_nlopt <- nloptr(x0 = c(mu_r, sigma_r), eval_f = logvraissemblance_pour_nloptr, opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 10000, tol_rel=1e-15, xtol_abs=1e-15), lb = c(1, 1))
    
    ## Mises à jours des paramètres
    mu_r <- res_nlopt$solution[1]
    sigma_r <- res_nlopt$solution[2]
    pi1_r <- T1/n
    pi2_r <- T2/n
    pi3_r <- T3/n
    
    ##Sauvegarde des historiques des paramètres d'interes
    suite_lambda <- c(suite_lambda, lambda_r)
    suite_mu <- c(suite_mu, mu_r)
    suite_sigma <- c(suite_sigma, sigma_r)
    suite_pi1 <- c(suite_pi1, pi1_r)
    suite_pi2 <- c(suite_pi2, pi2_r)
    suite_pi3 <- c(suite_pi3, pi3_r)
    
    
    if (length(Lvc)>2 & abs((tail(Lvc, 1) - tail(Lvc, 2)[1])) < 0.1) { break }
  })
  
  
  ## On vérifie qu'il n'y a pas eu une seule erreur de nlopt dans le repeat
  if(!(inherits(t, "try-error"))){
    
    ## On sauvegarde alors tous les paramètres d'intérêts finaux de notre algorithme EM
    Lvc <- c(Lvc, logvraissemblance(c(mu_r, sigma_r))) #Pour avoir autant de points que les paramètres
    ma_liste <- list(logvraisemblance  = Lvc, lambda = suite_lambda, mu = suite_mu, sigma = suite_sigma, pi1 = suite_pi1, pi2 = suite_pi2, pi3 = suite_pi3)
  }
  return (ma_liste)
}
toc()
#>Application: 1672.61 sec elapsed
stopCluster(nbre_coeurs_voulu)
save(resultats, file = "resultats_application.RData")
################################################################################
##--------------------------------RESULTATS-----------------------------------##
################################################################################

resultats <- load(file = "resultats_application.RData")

## On veut maintenant récupérer le meilleur des random starts, celui qui a la log-vraisemblance complétée la plus élevée.
max_index <- 1
max_lvc <- tail(resultats[[1]][[1]], 1)
for (i in 9:nbre_random_start){
  print(i)
  if (tail(resultats[[i]][[1]], 1)>max_lvc) {
    print("b")
    max_lvc <- tail(resultats[[i]][[1]], 1)
    max_index <- i
  }
}
max_index
max_lvc
max(resultats[[max_index]][[1]])
## Valeur de Theta
c(tail(resultats[[max_index]][[1]], 1), tail(resultats[[max_index]][[2]], 1), tail(resultats[[max_index]][[3]], 1), tail(resultats[[max_index]][[4]], 1), tail(resultats[[max_index]][[5]], 1), tail(resultats[[max_index]][[6]], 1), tail(resultats[[max_index]][[7]], 1))

