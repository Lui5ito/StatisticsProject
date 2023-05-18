library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(plotly)
library(cowplot)
library(DPQ)
library(rlist)


rm(list=ls())
set.seed(1664)
#--# Pour créer l'échantillon #--#
#####

load("gene_expression.RData")

sum_gene_df <- as.data.frame(colSums(gene_expression)) %>% 
  rename("nb_transcrit" =`colSums(gene_expression)`) %>% 
  filter(nb_transcrit != 0)

sum_gene_df$classe <- case_when(sum_gene_df$nb_transcrit < 2000 ~ "1",
                                sum_gene_df$nb_transcrit >= 2000 & sum_gene_df$nb_transcrit < 12000 ~ "2",
                                sum_gene_df$nb_transcrit >= 12000 ~ "3")

data <- sum_gene_df[,1]


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
    
    ## Le calcul de lambda_r est simple
    suite_lambda <- c(suite_lambda, sum(t1*data)/T1)
    
    
    ## Le calcul de mu_r est plus complexe, on a besoin d'intermediaires, et d'une fonction pour laquelle on va chercher le zéro.
    ## Le calcul de sigma_r est plus complexe, on a besoin d'intermediaires, et d'une fonction pour laquelle on va chercher le zéro.
    res <- optim(par = c(mu_r, sigma_r), fn = logvraissemblance, method = "L-BFGS-B", control = list(fnscale=-1), lower = c(1, 1))
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
    
    if (length(Lvc)>2 & abs((tail(Lvc, 1) - tail(Lvc, 2)[1])) < 0.1) { break }
  })
  
  if(!(inherits(t, "try-error"))){
    ma_liste <- list(logvraisemblance  = Lvc, lambda = suite_lambda, mu = suite_mu, sigma = suite_sigma, pi1 = suite_pi1, pi2 = suite_pi2, pi3 = suite_pi3)
    results <- list.append(results, ma_liste)
    print(i)
  }
}



plot(Lvc)
plot(suite_lambda)
plot(suite_mu)
plot(suite_sigma)
plot(suite_pi1)
plot(suite_pi2)
plot(suite_pi3)
c(tail(suite_lambda, 1), tail(suite_mu, 1), tail(suite_sigma, 1))
