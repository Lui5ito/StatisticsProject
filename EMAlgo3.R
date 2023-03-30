library(Rmpfr)
library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(plotly)
library(cowplot)

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

sum_gene_df %>%
  group_by(classe) %>%
  summarize(total = n())

echant <- rbind((sum_gene_df %>% filter(classe == 1))[sample(nrow(sum_gene_df %>% filter(classe == 1)), 500), ], 
                (sum_gene_df %>% filter(classe == 2))[sample(nrow(sum_gene_df %>% filter(classe == 2)), 500), ], 
                (sum_gene_df %>% filter(classe == 3))[sample(nrow(sum_gene_df %>% filter(classe == 3)), 500), ])


#####

data <- mpfr(echant[,1], 128)
data <- mpfr(sum_gene_df[,1], 128)
n <- length(data)

lambda_r <- 1
mu_r <- 5000
sigma_r <- 10000
pi1_r <- 1/3
pi2_r <- 1/3
pi3_r <- 1/3
Lvc <- NULL


repeat{
  phi1 <- dpois(data, lambda_r, log = FALSE)
  phi2 <- pnorm(data+1/2, mean = mu_r, sd = sigma_r) - pnorm(data-1/2, mean = mu_r, sd = sigma_r)
  phi3 <- pnorm(data+1/2, mean = 2*mu_r, sd = sqrt(2)*sigma_r) - pnorm(data-1/2, mean = 2*mu_r, sd = sqrt(2)*sigma_r)
  
  ln1 <- log(pi1_r*phi1)
  ln2 <- log(pi2_r*phi2)
  ln3 <- log(pi3_r*phi3)
  
  somme_phi_pondere <- pi1_r*phi1 + pi2_r*phi2 + pi3_r*phi3
  
  t1 <- pi1_r*phi1 / somme_phi_pondere
  t2 <- pi2_r*phi2 / somme_phi_pondere
  t3 <- pi3_r*phi3 / somme_phi_pondere
  
  lv1 <- sum(t1*ln1)
  lv2 <- sum(t2*ln2)
  lv3 <- sum(t3*ln3)
  
  T1 <- sum(t1)
  T2 <- sum(t2)
  T3 <- sum(t3)
  
  #--# On peut calculer la log-vraissemblance complétée à ce moment #--#
  LV_r <- lv1+lv2+lv3
  Lvc <- c(Lvc, as.numeric(LV_r))
  
  #--# On peut maintenant calculer theta_r #--#
  
  ## Le calcul de lambda_r est simple
  lambda_r <- sum(t1*data)/T1
  
  ## Le calcul de mu_r est plus complexe, on a besoin d'intermediaires, et d'une fonction pour laquelle on va chercher le zéro.
  dLv_dmu <- function(mu){
    gamma <- -dnorm(data+0.5, mean = mu, sd = sigma_r) + dnorm(data-0.5, mean = mu, sd = sigma_r)
    delta <- pnorm(data+0.5, mean = mu, sd = sigma_r) - pnorm(data-0.5, mean = mu, sd = sigma_r)
    khi <- gamma/delta
    psi <- sum(t2*khi)
    
    kappa <- -dnorm(data+0.5, mean = 2*mu, sd = sqrt(2)*sigma_r) + dnorm(data-0.5, mean = 2*mu, sd = sqrt(2)*sigma_r)
    tau <- pnorm(data+0.5, mean = 2*mu, sd = sqrt(2)*sigma_r) - pnorm(data-0.5, mean = 2*mu, sd = sqrt(2)*sigma_r)
    upsilon <- kappa/tau
    rho <- sum(t3*upsilon)
    
    resultat <- as.numeric(psi+rho)
    
    return (resultat)
  }
  mu_r <- uniroot(dLv_dmu, interval = c(10, 100000))$root
  
  ## Le calcul de sigma_r est plus complexe, on a besoin d'intermediaires, et d'une fonction pour laquelle on va chercher le zéro.
  dLv_dsigma <- function(sigma){
    gamma <- dnorm(data+0.5, mean = mu_r, sd = sigma)*((mu_r-(data+0.5))/sigma) - dnorm(data-0.5, mean = mu_r, sd = sigma)*((mu_r-(data-0.5))/sigma)
    delta <- pnorm(data+0.5, mean = mu_r, sd = sigma) - pnorm(data-0.5, mean = mu_r, sd = sigma)
    khi <- gamma/delta
    psi <- sum(t2*khi)
    
    kappa <- dnorm(data+0.5, mean = 2*mu_r, sd = sqrt(2)*sigma)*((2*mu_r-(data+0.5))/sqrt(2)*sigma) - dnorm(data-0.5, mean = 2*mu_r, sd = sqrt(2)*sigma)*((2*mu_r-(data-0.5))/sqrt(2)*sigma)
    tau <- pnorm(data+0.5, mean = 2*mu_r, sd = sqrt(2)*sigma) - pnorm(data-0.5, mean = 2*mu_r, sd = sqrt(2)*sigma)
    upsilon <- kappa/tau
    rho <- sum(t3*upsilon)
    
    resultat <- as.numeric(psi+rho)
    
    return (resultat)
  }
  sigma_r <- uniroot(dLv_dsigma, c(10, 100000))$root
  
  ##Les calculs des pi sont simples:
  pi1_r <- T1/n
  pi2_r <- T2/n
  pi3_r <- T3/n
  
  print(length(Lvc))
  
  if (length(Lvc)>2 & abs((tail(Lvc, 1) - tail(Lvc, 2)[1])) < 1) { break }
}































































