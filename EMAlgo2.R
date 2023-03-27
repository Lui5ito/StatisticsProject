###### Mise en place de l'algorithme EM avec Optim#####

EM <- function(data, lambda0, mu0, sigma0, pi1, pi2, pi3) {
  
  #On définit les variables qui vont varier [r]
  lambda_r <- lambda0
  mu_r <- mu0
  sigma_r <- sigma0
  pi1_r <- pi1
  pi2_r <- pi2
  pi3_r <- pi3
  
  #Variable useful throughout the algorithm
  n <- length(data)
  LVC <- NULL
  j <- 0
  
  ### Début du 'vrai' algorithme ### 
  
  repeat {
    t1 <- NULL
    t2 <- NULL
    t3 <- NULL
    ln1 <- NULL
    ln2 <- NULL
    ln3 <- NULL
    #On commence par calculer les tik(theta[r-1])

    phi1 <- function(lambda){
      return(dpois(data, lambda, log = FALSE))
    }
    phi2 <- function(mu, sigma){
      return(pnorm(data+1/2, mean = mu, sd = sigma) - pnorm(data-1/2, mean = mu, sd = sigma))
    }
    phi3 <- function(mu, sigma){
      return(pnorm(data+1/2, mean = 2*mu, sd = sqrt(2)*sigma) - pnorm(data-1/2, mean = 2*mu, sd = sqrt(2)*sigma))
    }
    
    somme_phi_pondere <- pi1_r*phi1(lambda_r) + pi2_r*phi2(mu_r, sigma_r) + pi3_r*phi3(2*mu_r, sqrt(2)*sigma_r)
    
    ln1 <- function(lambda){
      return(log(pi1_r)+log(phi1(lambda)))
    }
    
    ln2 <- function(mu, sigma){
      return(log(pi2_r)+log(phi2(mu, sigma)))
    }
    
    ln3 <- function(mu, sigma){
      return(log(pi3_r)+log(phi3(mu, sigma)))
    }
    
    t1 <- pi1_r*phi1(lambda_r) / somme_phi_pondere
    t2 <- pi2_r*phi2(mu_r, sigma_r) / somme_phi_pondere
    t3 <- pi3_r*phi3(mu_r, sigma_r) / somme_phi_pondere
    
    #On peut calculer la LVC de l'étape précédente
    LVC <- c(LVC, (sum(t1*ln1(lambda_r)) + sum(t2*ln2(mu_r, sigma_r)) + sum(t3*ln3(mu_r, sigma_r))))
    
    #On calculs les Tk
    T1 <- sum(t1)
    T2 <- sum(t2)
    T3 <- sum(t3)
    
    #On peut calculer les nouveaux PIk
    pi1_r <- T1/n
    pi2_r <- T2/n
    pi3_r <- T3/n
    
    #On calculs les paramètres des distributions (Etape M), ici on la fait avec optim
    to_argmax <- function(par){
      lambda <- exp(par[1])
      mu <- exp(par[2])
      sigma <- exp(par[3])
      a <- sum(t1*ln1(lambda)) + sum(t2*ln2(mu, sigma)) + sum(t3*ln3(mu, sigma))
      return(a)
    }
    
    print(c(lambda_r, mu_r, sigma_r))
    print(phi1(lambda_r))
    chapeau <- optim(par = c(lambda_r, mu_r, sigma_r), fn = to_argmax)
    
    lambda_r <- chapeau$par[1]
    mu_r <- chapeau$par[2]
    sigma_r <- chapeau$par[3]
    
    
    
    
    
    #Tant que la diff entre les deux dernier éléments de LVC est supérieur à 10^-3
    print(j)
    j <- j+1
    #if (abs((tail(LVC, 1) - tail(LVC, 2)[1])) < 1) {
    # break
    #}
    if (j>100) { break }
  }
  
  return (list(lambda = lambda_r, mu1 = mu1_r, sigma1 = sigma1_r, mu2 = mu2_r, sigma2 = sigma2_r, pi1 = pi1_r, pi2 = pi2_r, pi3 = pi3_r, LVC = LVC))
}


##### Test de mon algo EM #####
data <- c(rpois(100, 3), trunc(rnorm(100, mean = 5000, sd = 1.5)), trunc(rnorm(100, mean = 10000, sd = 3)))
data <- echant[,1]
test <- EM(data, 3, 500, 1, 1/3, 1/3, 1/3)
dpois(data, 3)
