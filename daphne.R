#PARAMETRES
library(CASdatasets)
data("freMPL10")
data = freMPL10
claimNB = data$ClaimNbResp+data$ClaimNbNonResp+data$ClaimNbParking+data$ClaimNbFireTheft+data$ClaimNbWindscreen
#claimNB=data$ClaimNbWindscreen+data$ClaimNbParking+data$ClaimNbFireTheft
#variable à expliquer
y = claimNB
#variables explicatives
x = cbind(data$BonusMalus, data$VehUsage, data$SocioCateg, data$RiskArea, data$Age,data$VehAge)
#contrainte
data$Age_num = as.numeric(data$Age)
z=data$Age_num


#Contrainte : HGR
library(acepack)
HGR_constraint <- function(beta, y, x, z) {
  hgr0 = ace(y,z)
  correlation <- as.numeric(cor(hgr0$tx, hgr0$ty,method='pearson' )^2)
  return(correlation)
}

HGR_optim <- function(beta){
  y = claimNB
  x = cbind(data$BonusMalus, data$VehUsage, data$SocioCateg, data$RiskArea, data$Age,data$VehAge)
  z = data$Age_num
  return(as.numeric(HGR_constraint(beta, y, x, z)))
}


#Fonction objectif avec contrainte
model_objective <- function(beta, y, x, z,lambda) {
  # Calcul de la fonction de vraisemblance
  L <- -sum(exp(x%*%beta))+sum(x%*%beta*y)-sum(log(factorial(y)))
  # Calcul de la contrainte de corrélation
  C <- -lambda*HGR_constraint(beta, y, x, z)
  # Fonction objectif totale
  return(L+C)
}

model_objective_nlopt <- function(beta){
  y = claimNB
  x = cbind(data$BonusMalus, data$VehUsage, data$SocioCateg, data$RiskArea, data$Age,data$VehAge)
  z = data$Age_num
  lambda = 1
  return(-model_objective(beta, y, x, z, lambda))
}

vraisemblance <- function(params){
  
}




library(nloptr)
# Estimation des paramètres
beta_init = runif(ncol(x), min = 0, max = 2)
model_objective_nlopt(beta_init)
res_nlopt <- nloptr(x0 = beta_init, eval_f = model_objective_nlopt, opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 10000, tol_rel=1e-15, xtol_abs=1e-15))
res_nlopt$solution
model_objective_nlopt(res_nlopt$solution)
beta_opt = result$par
lambda_opt = result$value - sum((y - x %*% beta_opt)^2)










