library(dplyr) #nécessaire pour l'importation des données
library(DPQ) #logspace.sub
library(nloptr) #nolptr
library(doParallel) #créer les clusters
library(foreach) #permet le parallélisme
library(tictoc)
library(ggplot2)

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

nbre_random_start <- 300 #### On prend plus que 2x le nombre de paramètres car on a un échantillon dispoportionné
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

resultats <- readRDS("resultats_application.RData")

## On veut maintenant récupérer le meilleur des random starts, celui qui a la log-vraisemblance complétée la plus élevée.
max_index <- 1
length(resultats[[8]])
max_lvc <- tail(resultats[[1]][[1]], 1)
for (i in 2:length(resultats)){
  if (length(resultats[[i]]) != 0){
    if (tail(resultats[[i]][[1]], 1)>max_lvc) {
      max_lvc <- tail(resultats[[i]][[1]], 1)
      max_index <- i
    }
  }
}

## Les listes se récupère avec 
resultats[[max_index]][[1]] ## pour la vraisemblance
ggplot(as.data.frame(resultats[[max_index]][[1]]), aes(x=seq_along(resultats[[max_index]][[1]]), y=resultats[[max_index]][[1]])) +
  geom_point(color = "#66CCFF", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab("log-vraisemblance complétée")
resultats[[max_index]][[2]] ## pour lambdahttp://127.0.0.1:31587/graphics/115581a2-66b6-4da8-b610-c06496e6b790.png
ggplot(as.data.frame(resultats[[max_index]][[2]]), aes(x=seq_along(resultats[[max_index]][[2]]), y=resultats[[max_index]][[2]])) +
  geom_point(color = "#66CCFF", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab(expression(lambda))
resultats[[max_index]][[3]] ## pour mu
ggplot(as.data.frame(resultats[[max_index]][[3]]), aes(x=seq_along(resultats[[max_index]][[3]]), y=resultats[[max_index]][[3]])) +
  geom_point(color = "#66CCFF", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab(expression(mu))
resultats[[max_index]][[4]] ## pour sigma
ggplot(as.data.frame(resultats[[max_index]][[4]]), aes(x=seq_along(resultats[[max_index]][[4]]), y=resultats[[max_index]][[4]])) +
  geom_point(color = "#66CCFF", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab(expression(sigma))
resultats[[max_index]][[5]] ## pour pi1
ggplot(as.data.frame(resultats[[max_index]][[5]]), aes(x=seq_along(resultats[[max_index]][[5]]), y=resultats[[max_index]][[5]])) +
  geom_point(color = "#66CCFF", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab(expression(pi[1]))
resultats[[max_index]][[6]] ## pour pi2
ggplot(as.data.frame(resultats[[max_index]][[6]]), aes(x=seq_along(resultats[[max_index]][[6]]), y=resultats[[max_index]][[6]])) +
  geom_point(color = "#66CCFF", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab(expression(pi[2]))
resultats[[max_index]][[7]] ## pour pi3
ggplot(as.data.frame(resultats[[max_index]][[7]]), aes(x=seq_along(resultats[[max_index]][[7]]), y=resultats[[max_index]][[7]])) +
  geom_point(color = "#66CCFF", size = 2.5) +
  theme_bw() +
  xlab("Itération") +
  ylab(expression(pi[3]))



## On récupère le theta_hat final pour pouvoir classer les individus
lambda_hat <- tail(resultats[[max_index]][[2]], 1)
mu_hat <- tail(resultats[[max_index]][[3]], 1)
sigma_hat <- tail(resultats[[max_index]][[4]], 1)
pi1_hat <- tail(resultats[[max_index]][[5]], 1)
pi2_hat <- tail(resultats[[max_index]][[6]], 1)
pi3_hat <- tail(resultats[[max_index]][[7]], 1)


## On peut générer un échantillon de notre loi de mélange et tracer son histogramme
taille_echantillon <- 350000
echantillon_hat <- c(rpois(taille_echantillon*pi1_hat, lambda_hat+3), round(rnorm(taille_echantillon*pi2_hat, mu_hat, sigma_hat)), round(rnorm(taille_echantillon*pi3_hat, 2*mu_hat, sqrt(2)*sigma_hat)))
echantillon_hat <- as.data.frame(cbind(echantillon_hat, c(rep(1, taille_echantillon*pi1_hat), rep(2, taille_echantillon*pi2_hat), rep(3, taille_echantillon*pi3_hat))))
colnames(echantillon_hat) <- c("echantillon", "groupe")
echantillon_hat <- echantillon_hat[which(echantillon_hat$echantillon > 0),]
echantillon_hat <- echantillon_hat[sample(nrow(echantillon_hat)),]

ggplot(echantillon_hat, aes(x = seq_along(echantillon), y = echantillon, color = as.factor(groupe))) +
  geom_point(shape = '.', alpha = 0.4) +
  scale_y_continuous(trans = 'log10') +
  theme_bw() +
  ggtitle("Scatter plot,\n and a first determination of the groups") +
  xlab("Index") + ylab("Number of transcript") +
  scale_color_manual(values=c("black", "#66CCFF", "red"))

ggplot(sum_gene_df, aes(x = seq_along(nb_transcrit), y = nb_transcrit)) +
  geom_point(shape = '.') +
  scale_y_continuous(trans = 'log10') +
  theme_bw() +
  ggtitle("Scatter plot,\n and a first determination of the groups") +
  xlab("Index") + ylab("Number of transcript")

## Groupe 1
echantillon_hat_1 <- echantillon_hat[which(echantillon_hat$groupe == 1),]
echantillon_hat_1 <- subset(echantillon_hat_1, select = -groupe)

data_classee_1 <- data_classee[which(data_classee$groupe == 1),]
data_classee_1 <- subset(data_classee_1, select = -groupe)

DF1 <- data.frame(variable=rep(c('Echantillon issu de la loi de mélange', 'Données réelles'), each=300000), value=c(echantillon_hat_1$echantillon[1:300000] , data_classee_1$data[1:300000]))

ggplot(DF1) + 
  stat_ecdf(aes(value, color=variable), size = 1) +
  scale_color_manual(values=c("red", "#66CCFF")) +
  theme_bw()

tapply(DF1$value, DF1$variable, summary)


## Groupe 2
echantillon_hat_2 <- echantillon_hat[which(echantillon_hat$groupe == 2),]
echantillon_hat_2 <- subset(echantillon_hat_2, select = -groupe)

data_classee_2 <- data_classee[which(data_classee$groupe == 2),]
data_classee_2 <- subset(data_classee_2, select = -groupe)

DF2 <- data.frame(variable=rep(c('Echantillon issu de la loi de mélange', 'Données réelles'), each=300000), value=c(echantillon_hat_2$echantillon[1:300000] , data_classee_2$data[1:300000]))

ggplot(DF2) + 
  stat_ecdf(aes(value, color=variable), size = 1) +
  scale_color_manual(values=c("red", "#66CCFF")) +
  theme_bw()

tapply(DF2$value, DF2$variable, summary)

## Groupe 3
echantillon_hat_3 <- echantillon_hat[which(echantillon_hat$groupe == 3),]
echantillon_hat_3 <- subset(echantillon_hat_3, select = -groupe)

data_classee_3 <- data_classee[which(data_classee$groupe == 3),]
data_classee_3 <- subset(data_classee_3, select = -groupe)

DF3 <- data.frame(variable=rep(c('Echantillon issu de la loi de mélange', 'Données réelles'), each=300000), value=c(echantillon_hat_3$echantillon[1:300000] , data_classee_3$data[1:300000]))

ggplot(DF3) + 
  stat_ecdf(aes(value, color=variable), size = 1) +
  scale_color_manual(values=c("red", "#66CCFF")) +
  theme_bw()

tapply(DF3$value, DF3$variable, summary)

###

DF <- data.frame(variable=rep(c('Echantillon issu de la loi de mélange', 'Données réelles'), each=length(echantillon_hat)), value=c(echantillon_hat$echantillon , data[1:length(echantillon_hat)]))

ggplot(DF) + 
  stat_ecdf(aes(value, color=variable), size = 1) +
  scale_color_manual(values=c("red", "#66CCFF")) +
  theme_bw()



  ggplot(as.data.frame(echantillon_hat), aes(x = echantillon_hat)) + 
  stat_ecdf(geom = "point",  color = "red") +
  labs(title="Cumulative Density Function \nfor Poisson",
       y = "F(x)", x="x") +
  theme_bw()

ggplot(as.data.frame(data), aes(x = data)) + 
  stat_ecdf(geom = "point",  color = "red") +
  labs(title="Cumulative Density Function \nfor Poisson",
       y = "F(x)", x="x") +
  theme_bw()

## Les probabilités sa caculent par le k qui maximise le tik(theta)
## Il faut donc calculer les tik(theta) pour les individus et tous les groupes
logphi1 <- dpois(data, lambda_hat, log = TRUE)
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

## On compte le nombre d'individu par groupe, et leur proportion (on retrouve les proportions du modeles)
population_par_classe <- count(data_classee, groupe)
population_par_classe$proportion <- round((population_par_classe$n / length(data))*100, 3)
population_par_classe

## On peut compter le nombre d'individu par groupe et leur proportion avant le modèle
pop_classe_avant <- count(sum_gene_df, classe)
pop_classe_avant$proportion <- round((pop_classe_avant$n / length(data))*100, 3)
pop_classe_avant


## On compte le nombre d'occurence de chaque individus et leur proportion par rapport à la population totale
## Piste pour amélioration, est ce que on peut retirer les premières gouttelettes ?
population_par_sumgene <- count(data_classee, data)
population_par_sumgene$proportion <- round((population_par_sumgene$n / length(data))*100, 3)
population_par_sumgene

plot_classe <- ggplot(data = data_classee) +
  geom_boxplot(aes(x = factor(groupe), y = data)) +
  scale_x_discrete(breaks = 1:3, labels = c("Groupe 1", "Groupe 2", "Groupe 3")) +
  #scale_y_continuous (trans='log10') +
  labs(title = "Etude de la performance des algorithmes d'optimisation",
       x = "Taille de l'échantillon",
       y = "Temps de convergence (échelle logarithmique)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  theme_bw()
plot_classe


