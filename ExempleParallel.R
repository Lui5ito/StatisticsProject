library(tictoc)
library(foreach)
library(doParallel)

rm(list=ls())
#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(0.5*cores[1]) #not to overload your computer
registerDoParallel(cl)

tic("no_parallel")

for (subject in 1:400){
  for (number in 1:5000000){
    sqrt(number)}
}

toc()
#> no_parallel: 271.312 sec elapsed

tic("parallel")

foreach (subject=1:400) %dopar% {
  for (number in 1:5000000){
    sqrt(number)}
}

toc()
#> parallel: 65.654 sec elapsed
#> 
#> #stop cluster
stopCluster(cl)
