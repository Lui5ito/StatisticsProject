library(DPQ)
library(nloptr)
library(doParallel)
library(foreach)

#########################################################
#####DONNE LE BON RESULTAT AVEC LES DONNEES SIMULEES#####
#########################################################

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
##-------------------------------ALGORITHME-----------------------------------##
################################################################################

## Initialisation du parallélisme
nbre_coeurs_disponibles = detectCores()
pourcentage_des_coeurs_voulu <- 0.5
nbre_coeurs_voulu <- makeCluster(pourcentage_des_coeurs_voulu*nbre_coeurs_disponibles)
registerDoParallel(nbre_coeurs_voulu)

## Initialisation des paramètres des lois cible
lambda_cible <- 3
mu_cible <- 5050
sigma_cible <- 1010


taille_echantillon <- c(10^2)
nbre_repetition <- 2
algorithmes <- c("NLOPT_LN_NELDERMEAD", "NLOPT_LN_COBYLA", "NLOPT_LN_BOBYQA") #### Pas de "NLOPT_LN_SBPLX" car trop proche de L-BFGS-B et difficile de justifier qu'il ne marche pas.  
nbre_random_start <- 4 #### 2x le nombre de paramètres.
nbre_parametre_interet <- 9 #### log-vraisemblance complétée, temps, iterations, lambda, mu, sigma, pi1, pi2, pi3

#!!# On ne s'intéresse ici que aux resultats finaux. On ne s'intéresse pas à l'évolution de la log-vraisemblance complétée ni à l'évolution de theta.
#!!# En revanche, ces aspects sont fondamentals dans le rapport et dans l'exécution de l'algorithme sur nos données réelles.

resultats <- array(data = NA, c(nbre_parametre_interet, length(taille_echantillon), nbre_repetition, length(algorithmes), nbre_random_start))
# resultats[1:9,i,j,k,l]

for (i in 1:length(taille_echantillon)){
  taille_echantillon_en_cours <- taille_echantillon[i]
  
  for (j in 1:nbre_repetition) {
    data <- c(rpois(taille_echantillon_en_cours/2, lambda_cible), round(rnorm(taille_echantillon_en_cours/4, mu_cible, sigma_cible)), round(rnorm(taille_echantillon_en_cours/4, 2*mu_cible, sqrt(2)*sigma_cible)))
    
    for (k in 1:length(algorithmes)){
      algo_en_cours <- algorithmes[k]
      
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
          resultats[1, i, j, k, l] <- logvraissemblance(c(mu_r, sigma_r)) ## log-vraisemblance complétée
          resultats[2, i, j, k, l] <- mean(temps_nlopt) ## On prend la moyenne de tous les nlopt utilisés lors de cet algorithme EM
          resultats[3, i, j, k, l] <- mean(iteration_nlopt) ## On doit encore prendre la moyenne du nombre d'itération de nlopt lors de l'algorithme EM
          resultats[4, i, j, k, l] <- lambda_r
          resultats[5, i, j, k, l] <- mu_r
          resultats[6, i, j, k, l] <- sigma_r
          resultats[7, i, j, k, l] <- pi1_r
          resultats[8, i, j, k, l] <- pi2_r
          resultats[9, i, j, k, l] <- pi3_r
        }
        
      }
    }
  }
}



summary(resultats)

## On récupère un tableau juste avec la log-vraisemblance complétée
ll <- resultats[1, , , , ]

## On récupère le meilleur des random-start : 
apply(X = resultats, MARGIN = c(1, 2, 3, 4), FUN = max, na.rm=TRUE)[1,]


## Moyenne de log-vraisemblance complétée des trois algorithmes
apply(X = resultats, MARGIN = c(1,4), FUN = mean, na.rm=TRUE)[1,]

## Maximum de log-vraisemblance complétée des trois algorithmes
apply(X = resultats, MARGIN = c(1,4), FUN = max, na.rm=TRUE)[1,]

## Moyenne de lambda des trois algorithmes
apply(X = resultats, MARGIN = c(1,4), FUN = mean, na.rm=TRUE)[4,]































for (k in nbre_donnes){
  data <- c(rpois(k/2, lambda_cible), round(rnorm(k/4, mu_cible, sigma_cible)), round(rnorm(k/4, 2*mu_cible, sqrt(2)*sigma_cible)))
  
  for (j in algos) {
    ma_liste <- list()
    results <- list()
    
    for (i in 1:10) {
      
      suite_lambda <- NULL
      suite_mu <- NULL
      suite_sigma <- NULL
      suite_pi1 <- NULL
      suite_pi3 <- NULL
      suite_pi2 <- NULL
      ma_liste <- NULL
      
      n <- length(data)
      lambda_r <- data[sample(1:n, 1)]
      mu_r <- data[sample(1:n, 1)]
      sigma_r <- 1000*sd(data)
      pi1_r <- 1/3
      pi2_r <- 1/3
      pi3_r <- 1/3
      Lvc <- NULL
      suite_lambda <- c(lambda_r)
      suite_sigma <- c(sigma_r)
      suite_mu <- c(mu_r)
      suite_pi1 <- c(pi1_r)
      suite_pi2 <- c(pi2_r)
      suite_pi3 <- c(pi3_r)
      
      t <- try(repeat{
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
        
        #--# On peut calculer la log-vraissemblance complétée à ce moment #--#
        LV_r <- lv1+lv2+lv3
        Lvc <- c(Lvc, LV_r)
        #--# On peut maintenant calculer theta_r #--#
        lambda_r <- sum(t1*data)/T1
        suite_lambda <- c(suite_lambda, sum(t1*data)/T1)
        
        #res <- nloptr(x0 = c(mu_r, sigma_r), eval_f = logvraissemblance_pour_nloptr, opts = list(algorithm = "NLOPT_LN_NELDERMEAD", maxeval = 10000, tol_rel=1e-15, xtol_abs=1e-15), lb = c(1, 1))
        
        res <- nloptr(x0 = c(mu_r, sigma_r), eval_f = logvraissemblance_pour_nloptr, opts = list(algorithm = j, maxeval = 10000, tol_rel=1e-15, xtol_abs=1e-15), lb = c(1, 1))
        
        
        ##Mises à jours des paramètres
        
        mu_r <- res$solution[1]
        sigma_r <- res$solution[2]
        pi1_r <- T1/n
        pi2_r <- T2/n
        pi3_r <- T3/n
        
        ##Sauvegarde des historiques des paramètres
        
        suite_mu <- c(suite_mu, res$solution[1])
        suite_sigma <- c(suite_sigma, res$solution[2])
        suite_pi1 <- c(suite_pi1, pi1_r)
        suite_pi2 <- c(suite_pi2, pi2_r)
        suite_pi3 <- c(suite_pi3, pi3_r)
        
        if (length(Lvc)>2 & abs((tail(Lvc, 1) - tail(Lvc, 2)[1])) < 0.1) { break }
      })
      
      if(!(inherits(t, "try-error"))){
        Lvc <- c(Lvc, logvraissemblance(c(mu_r, sigma_r))) #Pour avoir autant de points que les paramètres
        ma_liste <- list(logvraisemblance  = tail(Lvc, 1), iteration = res$iterations, lambda = tail(suite_lambda, 1), mu = tail(suite_mu, 1), sigma = tail(suite_sigma, 1), pi1 = tail(suite_pi1, 1), pi2 = tail(suite_pi2, 1), pi3 = tail(suite_pi3, 1))
        results <- list.append(results, ma_liste)
        erreurs_optim <- c(erreurs_optim, res$message)
      }
    }
    logvrais_final <- c()
    for (i in 1:length(results)){
      logvrais_final <- c(logvrais_final, max(results[[i]]$logvraisemblance))
    }
    index_max <- which.max(logvrais_final)
    ma_liste_algo <- list(algo = j, logvrais = tail(results[[index_max]]$logvraisemblance, 1), iteration = tail(results[[index_max]]$iteration, 1), lambda = tail(results[[index_max]]$lambda, 1), mu = tail(results[[index_max]]$mu, 1), sigma = tail(results[[index_max]]$sigma, 1), pi1 = tail(results[[index_max]]$pi1, 1), pi2 = tail(results[[index_max]]$pi2, 1), pi3 = tail(results[[index_max]]$pi3, 1))
    results_algos <- list.append(results_algos, ma_liste_algo)
    if (j == "NLOPT_LN_SBPLX"){iteration_SBPLX <- c(iteration_SBPLX, tail(results[[index_max]]$iteration, 1)) }
    if (j == "NLOPT_LN_NELDERMEAD"){iteration_NELDERMEAD <- c(iteration_NELDERMEAD, tail(results[[index_max]]$iteration, 1)) }
    if (j == "NLOPT_LN_COBYLA"){iteration_COBYLA <- c(iteration_COBYLA, tail(results[[index_max]]$iteration, 1)) }
    if (j == "NLOPT_LN_BOBYQA"){iteration_BOBYQA <- c(iteration_BOBYQA, tail(results[[index_max]]$iteration, 1)) }
  }
  
}
df <- as.data.frame(cbind(x = nbre_donnes, bobyqa = iteration_BOBYQA, cobyla = iteration_COBYLA, neldermead = iteration_NELDERMEAD, sbplx = iteration_SBPLX))

ggplot(df, aes(nbre_donnes)) +  
  geom_line(aes(y = bobyqa), color = "purple") +
  geom_point(aes(y = bobyqa), color = "purple") +
  geom_line(aes(y = cobyla), color = "red") +
  geom_point(aes(y = cobyla), color = "red") +
  geom_line(aes(y = neldermead), color = "green") +
  geom_point(aes(y = neldermead), color = "green") +
  geom_line(aes(y = sbplx), color = "blue") +
  geom_point(aes(y = sbplx), color = "blue")
