logvraissemblance <- function(lambda, mu, sigma, pi1, pi2, pi3){
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
  
  LV_r <- lv1+lv2+lv3
}

gradient_mu_sigma <- function(mu, sigma){
  return(c(dLv_dmu(mu), dLv_dsigma(sigma)))
}

dLv_dmu <- function(mu){
  gamma <- -dnorm(data+0.5, mean = mu, sd = sigma_r) + dnorm(data-0.5, mean = mu, sd = sigma_r)
  delta <- pnorm(data+0.5, mean = mu, sd = sigma_r) - pnorm(data-0.5, mean = mu, sd = sigma_r)
  khi <- gamma/delta
  psi <- sum(t2*khi)
  
  kappa <- -dnorm(data+0.5, mean = 2*mu, sd = sqrt(2)*sigma_r) + dnorm(data-0.5, mean = 2*mu, sd = sqrt(2)*sigma_r)
  tau <- pnorm(data+0.5, mean = 2*mu, sd = sqrt(2)*sigma_r) - pnorm(data-0.5, mean = 2*mu, sd = sqrt(2)*sigma_r)
  upsilon <- kappa/tau
  rho <- sum(t3*upsilon)
  
  return (psi+rho)
}

dLv_dsigma <- function(sigma){
  gamma <- dnorm(data+0.5, mean = mu_r, sd = sigma)*((mu_r-(data+0.5))/sigma) - dnorm(data-0.5, mean = mu_r, sd = sigma)*((mu_r-(data-0.5))/sigma)
  delta <- pnorm(data+0.5, mean = mu_r, sd = sigma) - pnorm(data-0.5, mean = mu_r, sd = sigma)
  khi <- gamma/delta
  psi <- sum(t2*khi)
  
  kappa <- dnorm(data+0.5, mean = 2*mu_r, sd = sqrt(2)*sigma)*((2*mu_r-(data+0.5))/sqrt(2)*sigma) - dnorm(data-0.5, mean = 2*mu_r, sd = sqrt(2)*sigma)*((2*mu_r-(data-0.5))/sqrt(2)*sigma)
  tau <- pnorm(data+0.5, mean = 2*mu_r, sd = sqrt(2)*sigma) - pnorm(data-0.5, mean = 2*mu_r, sd = sqrt(2)*sigma)
  upsilon <- kappa/tau
  rho <- sum(t3*upsilon)
  
  return (psi+rho)
}

pnorm_perso <- function(x){
  normal01 <- dnorm(x, mean = 5000, sd = 10000)
  p <- 0.2316419
  t <- 1/(1+p*x)
  t15 <- c(t, t**2, t**3, t**4, t**5)
  b <- c(0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429)
  c <- sum(t15*b)
  
  return(1+normal01*c)
}

