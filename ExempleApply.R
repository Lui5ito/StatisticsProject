ntaillech <- 3; nrep <- 10; nbalgo <- 4; nmulti <- 5
A <- array(runif(ntaillech*nrep*nbalgo*nmulti),c(ntaillech,nrep,nbalgo,nmulti))

B <- apply(A,1:3,max)

# Pour les moyennes et écart-type du tableau
C <- apply(B,c(1,3),mean)
C <- apply(B,c(1,3),sd)

# Pour les boxplots
library(reshape2)
D <- melt(B,varnames=c("ntailleech","nrep","nbalog"))

A[apply(A,1:3,which.max)]
