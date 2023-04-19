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
