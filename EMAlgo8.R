library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(plotly)
library(cowplot)
library(DPQ)
library(rlist)
library(nloptr)

#########################################################
#####DONNE LE BON RESULTAT AVEC LES DONNEES SIMULEES#####
#########################################################

rm(list=ls())
set.seed(1664)
#--# Pour créer l'échantillon #--#
#####

lambda_cible <- 4
mu_cible <- 7900
sigma_cible <- 800
data <- c(rpois(1000, lambda_cible), round(rnorm(500, mu_cible, sigma_cible)), round(rnorm(300, 2*mu_cible, sqrt(2)*sigma_cible)))

################################################################################
##----------------------------Test avec NLOPTR--------------------------------##
################################################################################



test <- nloptr(c(5000, 1000), logvraissemblance, opts = list(algorithm = "NLOPT_LN_NELDERMEAD", maxeval = 10000, tol_rel=1e-15, xtol_abs=1e-15), lb = c(1, 1), ub = c(10000, 10000))

nloptr()
test$status
test$message
test$solution





################################################################################
##-------------------------------FUNCTIONS------------------------------------##
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
results <- list()
erreurs_optim <- c()
for (i in 1:5) {
  
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
    
    ## Le calcul de lambda_r est simple
    suite_lambda <- c(suite_lambda, sum(t1*data)/T1)
    
    

    #res <- optim(par = c(mu_r, sigma_r), fn = logvraissemblance, method = "L-BFGS-B", control = list(fnscale=-1, ndeps=c(1e-8, 1e-8), pgtol = 1e-8, factr =  1e12), lower = c(1, 1))
    res <- nloptr(x0 = c(mu_r, sigma_r), eval_f = logvraissemblance_pour_nloptr, opts = list(algorithm = "NLOPT_LN_NELDERMEAD", maxeval = 10000, tol_rel=1e-15, xtol_abs=1e-15), lb = c(1, 1))
    suite_mu <- c(suite_mu, res$solution[1])
    suite_sigma <- c(suite_sigma, res$solution[2])
    
    ##Mises à jours des paramètres
    lambda_r <- sum(t1*data)/T1
    mu_r <- res$solution[1]
    sigma_r <- res$solution[2]
    pi1_r <- T1/n
    pi2_r <- T2/n
    pi3_r <- T3/n
    
    suite_pi1 <- c(suite_pi1, pi1_r)
    suite_pi2 <- c(suite_pi2, pi2_r)
    suite_pi3 <- c(suite_pi3, pi3_r)
    
    if (length(Lvc)>2 & abs((tail(Lvc, 1) - tail(Lvc, 2)[1])) < 0.1) { break }
  })
  
  if(!(inherits(t, "try-error"))){
    ma_liste <- list(logvraisemblance  = Lvc, lambda = suite_lambda, mu = suite_mu, sigma = suite_sigma, pi1 = suite_pi1, pi2 = suite_pi2, pi3 = suite_pi3)
    results <- list.append(results, ma_liste)
    erreurs_optim <- c(erreurs_optim, res$message)
    print(i)
  }
}


erreurs_optim


logvrais_final <- c()
for (i in 1:length(results)){
  logvrais_final <- c(logvrais_final, max(results[[i]]$logvraisemblance))
}
index_max <- which.max(logvrais_final)
plot(results[[index_max]]$logvraisemblance)
plot(results[[index_max]]$lambda)
plot(results[[index_max]]$mu)
plot(results[[index_max]]$sigma)
plot(results[[index_max]]$pi1)
plot(results[[index_max]]$pi2)
plot(results[[index_max]]$pi3)
c(tail(results[[index_max]]$lambda, 1), tail(results[[index_max]]$mu, 1), tail(results[[index_max]]$sigma, 1))





plot(Lvc)
plot(suite_lambda)
plot(suite_mu)
plot(suite_sigma)
plot(suite_pi1)
plot(suite_pi2)
plot(suite_pi3)
c(tail(suite_lambda, 1), tail(suite_mu, 1), tail(suite_sigma, 1))


list_mu <- seq(1:10000)
list_vrais_mu <- NULL
for (i in list_mu){
  list_vrais_mu <- c(list_vrais_mu, logvraissemblance_sigma1000(i))
}
plot(list_mu, list_vrais_mu)


logvraissemblance_mu5000 <- function(sigma) {
  return (logvraissemblance(c(5000, sigma)))
}
list_sigma <- seq(1:10000)
list_vrais_sigma <- NULL
for (i in list_mu){
  list_vrais_sigma <- c(list_vrais_sigma, logvraissemblance_mu5000(i))
}
plot(list_sigma, list_vrais_sigma)

data <- rnorm(100, mean = 0, sd = 1)
plot(data, pnorm(data, mean = 0, sd = 1))
