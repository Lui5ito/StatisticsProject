library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(plotly)
library(cowplot)
library(DPQ)

rm(list=ls())
set.seed(1664)
#--# Pour créer l'échantillon #--#
#####

lambda_cible <- 100
mu_cible <- 5000
sigma_cible <- 100
data <- c(rpois(700, lambda_cible), round(rnorm(200, mu_cible, sigma_cible)), round(rnorm(100, 2*mu_cible, sigma_cible)))

n <- length(data)
lambda_r <- 1 #data[sample(1:n, 1)]
mu_r <- 3000 #data[sample(1:n, 1)]
sigma_r <- 100*sd(data)
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


################################################################################
##-------------------------------ALGORITHME-----------------------------------##
################################################################################

repeat{
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
  print(Lvc)
  #--# On peut maintenant calculer theta_r #--#
  
  ## Le calcul de lambda_r est simple
  suite_lambda <- c(suite_lambda, sum(t1*data)/T1)
  
  
  ## Le calcul de mu_r est plus complexe, on a besoin d'intermediaires, et d'une fonction pour laquelle on va chercher le zéro.
  #mu_r <- uniroot(dLv_dmu, interval = c(10, 100000))$root
  #suite_mu <- c(suite_mu, mu_r)
  
  ## Le calcul de sigma_r est plus complexe, on a besoin d'intermediaires, et d'une fonction pour laquelle on va chercher le zéro.
  #sigma_r <- uniroot(dLv_dsigma, c(10, 100000))$root
  #suite_sigma <- c(suite_sigma, sigma_r)
  
  res <- optim(par = c(mu_r, sigma_r), fn = logvraissemblance, method = "L-BFGS-B", control = list(trace=6, fnscale=-1), lower = c(1, 1))
  suite_mu <- c(suite_mu, res$par[1])
  suite_sigma <- c(suite_sigma, res$par[2])
  
  ##Mises à jours des paramètres
  lambda_r <- sum(t1*data)/T1
  mu_r <- res$par[1]
  sigma_r <- res$par[2]
  pi1_r <- T1/n
  pi2_r <- T2/n
  pi3_r <- T3/n
  
  suite_pi1 <- c(suite_pi1, pi1_r)
  suite_pi2 <- c(suite_pi2, pi2_r)
  suite_pi3 <- c(suite_pi3, pi3_r)
  
  print(length(Lvc))
  
  if (length(Lvc)>2 & abs((tail(Lvc, 1) - tail(Lvc, 2)[1])) < 0.1) { break }
}

plot(Lvc)
plot(suite_lambda)
plot(suite_mu)
plot(suite_sigma)
plot(suite_pi1)
plot(suite_pi2)
plot(suite_pi3)
c(tail(suite_lambda, 1), tail(suite_mu, 1), tail(suite_sigma, 1))
