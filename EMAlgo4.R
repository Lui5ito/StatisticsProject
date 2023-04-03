library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(plotly)
library(cowplot)
library(DPQ)

source("Functions.R")

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

echant <- rbind((sum_gene_df %>% filter(classe == 1))[sample(nrow(sum_gene_df %>% filter(classe == 1)), 100), ], 
                (sum_gene_df %>% filter(classe == 2))[sample(nrow(sum_gene_df %>% filter(classe == 2)), 100), ], 
                (sum_gene_df %>% filter(classe == 3))[sample(nrow(sum_gene_df %>% filter(classe == 3)), 100), ])


#####

data <- echant[,1]
data <- sum_gene_df[,1]
n <- length(data)

lambda_r <- 1
mu_r <- 5000
sigma_r <- 1000
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
  
  #--# On peut maintenant calculer theta_r #--#
  
  ## Le calcul de lambda_r est simple
  suite_lambda <- c(suite_lambda, sum(t1*data)/T1)
  
  ## Le calcul de mu_r est plus complexe, on a besoin d'intermediaires, et d'une fonction pour laquelle on va chercher le zéro.
  #mu_r <- uniroot(dLv_dmu, interval = c(10, 100000))$root
  #suite_mu <- c(suite_mu, mu_r)
  
  ## Le calcul de sigma_r est plus complexe, on a besoin d'intermediaires, et d'une fonction pour laquelle on va chercher le zéro.
  #sigma_r <- uniroot(dLv_dsigma, c(10, 100000))$root
  #suite_sigma <- c(suite_sigma, sigma_r)
  
  optim(par = c(mu_r, sigma_r), fn = logvraissemblance, gr = gradient_mu_sigma, tail(suite_lambda, 1), tail(suite_pi1, 1), tail(suite_pi2, 1), tail(suite_pi3, 1), method = "L-BFGS-B", lower = c(1, 1), upper = c(20000, 20000))
  
  
  ##Les calculs des pi sont simples:
  pi1_r <- T1/n
  pi2_r <- T2/n
  pi3_r <- T3/n
  
  print(length(Lvc))
  
  if (length(Lvc)>2 & abs((tail(Lvc, 1) - tail(Lvc, 2)[1])) < 1) { break }
}

plot(Lvc)
plot(suite_lambda)
plot(suite_mu)
plot(suite_sigma)






























































