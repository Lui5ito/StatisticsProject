###### Mise en place de l'algorithme EM #####

EM <- function(data, lambda, mu1, sigma1, mu2, sigma2, pi1, pi2, pi3) {
  
  #On définit les variables qui vont varier [r]
  lambda_r <- lambda
  mu1_r <- mu1
  sigma1_r <- sigma1
  mu2_r <- mu2
  sigma2_r <- sigma2
  pi1_r <- pi1
  pi2_r <- pi2
  pi3_r <- pi3
  
  #Variable useful throughout the algorithm
  n <- length(data)
  t1 <- NULL
  t2 <- NULL
  t3 <- NULL
  ln1 <- NULL
  ln2 <- NULL
  ln3 <- NULL
  LVC <- c(0, 0.1)
  
  ### Début du 'vrai' algorithme ### 
  
  repeat {
    #On commence par calculer les tik(theta[r-1])
    for (i in 1:n) {
      phi1 <- ((lambda**data[i])*exp(-lambda))/factorial(data[i])
      phi2 <- (1/(data[i]*sqrt(2*sigma1)))*exp(-((log(data[i])-mu1)**2)/2*sigma1)
      phi3 <- (1/(data[i]*sqrt(2*pi*sigma2)))*exp(-((log(data[i])-mu2)**2)/2*sigma2)
      
      somme_phi_pondere <- pi1*phi1 + pi2*phi2 + pi3*phi3
      
      ln1 <- c(ln1, log(pi1*phi1))
      ln2 <- c(ln2, log(pi2*phi2))
      ln3 <- c(ln3, log(pi3*phi3))
      
      t1 <- c(t1, pi1*phi1 / somme_phi_pondere)
      t2 <- c(t2, pi2*phi2 / somme_phi_pondere)
      t3 <- c(t3, pi3*phi3 / somme_phi_pondere)
    }
    
    #On peut calculer la LVC de l'étape précédente
    LVC <- c(LVC, sum(t1*ln1) + sum(t2*ln2) + sum(t3*ln3))
    
    #On calculs les Tk
    T1 <- sum(t1)
    T2 <- sum(t2)
    T3 <- sum(t3)
    
    #On peut calculer les nouveaux PIk
    pi1_r <- T1/n
    pi2_r <- T2/n
    pi3_r <- T3/n
    
    #On calculs les paramètres des distributions
    lambda_r <- sum(t1*data)/T1
    mu1_r <- sum(t2*log(data))/T2
    temporaire <- (log(data)-mu1_r)**2
    sigma1_r <- sum(t2*temporaire)/T2
    mu2_r <- sum(t3*log(data))/T3
    temporaire <- (log(data)-mu2_r)**2
    sigma2_r <- sum(t3*temporaire)/T3
    
    #Tant que la diff entre les deux dernier éléments de LVC est supérieur à 10^-3
    if ((tail(LVC, 1) - tail(LVC, 2)[1]) < 0.1) {
      break
    }
  }
  
  return (list(lambda = lamnda_r, mu1 = mu1_r, sigma1 = sigma1_r, mu2 = mu2_r, sigma2 = sigma2_r, pi1 = pi1_r, pi2 = pi2_r, pi3 = pi3_r, LVC = LVC))
}


##### Test de mon algo EM #####
data <- c(rpois(100, 1), rlnorm(100, 0, 3/2), rlnorm(100, 2, 1))

test <- EM(data, 0.5, 1, 1, 1, 2, 1/3, 1/3, 1/3)
