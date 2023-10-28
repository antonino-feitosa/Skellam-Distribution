

library(boot)
library(skellam)
library(parallel)

start_time <- Sys.time()

set.seed(0)

lambda1 <- 4
lambda2 <- 5

lowerParametricSpace <- c(0.0001, 0.0001)
upperParametricSapce <- c(Inf, Inf)

n <- 100 # 25, 50, 75, 100
repetitions <- 100
repetitionsBoot <- 100

numcpus <- 3

jointDensity <- function(parans, sample){
  result <- -sum(dskellam(sample, parans[1], parans[2], log = TRUE))
  return (result)
}

momentEstimator <- function(sample){
  # https://www.emis.de/journals/BMMSS/pdf/acceptedpapers/2008-12-05a.pdf
  # https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwitgYTn8JiCAxX9pZUCHSspD3sQFnoECAsQAQ&url=http%3A%2F%2Fwww2.stat-athens.aueb.gr%2F~karlis%2FTR101_Soc2.pdf&usg=AOvVaw1btx41d_xPkIwoAa5iM-YJ&opi=89978449
  n <-length(sample) 
  v <- (n-1)/n * var(sample)
  m <- mean(sample)
  par1 <- (v + m)/2
  par2 <- (v - m)/2
  if(par1 < 0) par1 = 0
  if(par2 < 0) par2 = 0
  return (c(par1, par2))
}

maximumLikelihoodEstimator <- function(sample){
  result <- optim(c(1,1), jointDensity, gr = NULL, sample, method = "L-BFGS-B", lower=lowerParametricSpace, upper=upperParametricSapce)
  if(result$convergence != 0)
    result$par = c(NA,NA)
  result$par
}

boostrap <- function(sampledList, estimator){
  sampledList <- simplify2array(sampledList)
  parans <- estimator(sampledList)
  
  sampledListBoot <- replicate(repetitionsBoot, sample(simplify2array(sampledList), replace = TRUE), simplify = FALSE)
  result <- sapply(sampledListBoot, maximumLikelihoodEstimator)
  p1 <- 2 * parans[1] - mean(result[1,], na.rm = TRUE)
  p2 <- 2 * parans[2] - mean(result[2,], na.rm = TRUE)
  return (c(p1, p2))
}

resumeMonteCarlo <- function(name, sample, real){
  m <- mean(sample, na.rm = TRUE)
  v <- var(sample, na.rm = TRUE)
  b <- m - real
  e <- b^2 + v
  
  result <- paste("Estimator", name)
  result <- paste(result, "\tMean:", round(m,4) );
  result <- paste(result, "\tVar :", round(v,4) );
  result <- paste(result, "\tBias :", round(b,4) );
  result <- paste(result, "\tMSE :", round(e,4) );
  
  cat(result)
}

sampledList <- replicate(repetitions, rskellam(n, lambda1, lambda2), simplify = FALSE)
monteCarlo.mle <- simplify2array(mclapply(sampledList, maximumLikelihoodEstimator, mc.cores = numcpus))
monteCarlo.mme <- simplify2array(mclapply(sampledList, momentEstimator, mc.cores = numcpus))
monteCarlo.boot.mme <- simplify2array(mclapply(sampledList, boostrap, momentEstimator, mc.cores = numcpus))
monteCarlo.boot.mle <- simplify2array(mclapply(sampledList, boostrap, maximumLikelihoodEstimator, mc.cores = numcpus))

resumeMonteCarlo("mean(mle)     ", monteCarlo.mle[1,], lambda1)
resumeMonteCarlo("variance(mle) ", monteCarlo.mle[2,], lambda2)
resumeMonteCarlo("mean(mme)     ", monteCarlo.mme[1,], lambda1)
resumeMonteCarlo("variance(mme) ", monteCarlo.mme[2,], lambda2)
resumeMonteCarlo("mean(boot.mle)    ", monteCarlo.boot.mle[1,], lambda1)
resumeMonteCarlo("variance(boot.mle)", monteCarlo.boot.mle[2,], lambda2)
resumeMonteCarlo("mean(boot.mme)    ", monteCarlo.boot.mme[1,], lambda1)
resumeMonteCarlo("variance(boot.mme)", monteCarlo.boot.mme[2,], lambda2)


end_time <- Sys.time()
end_time - start_time


