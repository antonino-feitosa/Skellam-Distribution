
library(boot)
library(skellam)
library(parallel)


## https://en.wikipedia.org/wiki/Handball_at_the_Summer_Olympics
# 24.53 mean(c(21, 19, 23, 18, 32, 22, 27, 28, 26, 28, 22, 28, 25))
# 21.84 mean(c(16, 15, 22, 17, 25, 20, 26, 26, 24, 23, 21, 26, 23))

## https://en.wikipedia.org/wiki/Football_at_the_Summer_Olympics 1936+
# mean(c(2,3,2, 1, 3, 2, 4, 2, 3, 1,2,2, 3,3,2, 1,1,2, 1,2))
# mean(c(1, 1, 0, 0, 1, 1, 1,1,1, 0, 0, 1, 2,2,2, 0,0,1, 1,1 ))

## https://en.wikipedia.org/wiki/Basketball_at_the_Summer_Olympics
#mean(c(19,65,36, 89, 91, 73, 65, 51, 95, 86, 96, 76, 117, 95, 85, 84, 118, 107, 96, 87))
#mean(c(8, 21, 25, 55, 57, 59, 50, 50, 74, 77, 65, 63, 85, 69, 75, 69, 107, 100, 66, 82))

parans <- list(c(10,10), c(10,100), c(100,10), c(100,100))
sampleSizes <- c(25, 50, 75, 100)

experiments <- list(
  c(25, 10, 10),
  c(25, 100, 10),
  c(25, 10, 100),
  c(25, 100, 100),
  
  c(50, 10, 10),
  c(50, 100, 10),
  c(50, 10, 100),
  c(50, 100, 100),
  
  c(75, 10, 10),
  c(75, 100, 10),
  c(75, 10, 100),
  c(75, 100, 100),
  
  c(100, 10, 10),
  c(100, 100, 10),
  c(100, 10, 100),
  c(100, 100, 100)
)


start <- c(50,50)
lowerParametricSpace <- c(0.01, 0.01)
upperParametricSapce <- c(Inf, Inf)

repetitions <- 2000
repetitionsBoot <- 500

numcpus <- 3


jointDensity <- function(parans, sample){
  result <- -sum(dskellam(sample, parans[1], parans[2], log = TRUE))
  return (result)
}

momentEstimator <- function(sample){
  # https://www.emis.de/journals/BMMSS/pdf/acceptedpapers/2008-12-05a.pdf
  # https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwitgYTn8JiCAxX9pZUCHSspD3sQFnoECAsQAQ&url=http%3A%2F%2Fwww2.stat-athens.aueb.gr%2F~karlis%2FTR101_Soc2.pdf&usg=AOvVaw1btx41d_xPkIwoAa5iM-YJ&opi=89978449
  # https://uregina.ca/~kozdron/Teaching/Regina/252Winter16/Handouts/ch5.pdf (justificativa da var)
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
  result <- optim(start, jointDensity, gr = NULL, sample, method = "L-BFGS-B", lower=lowerParametricSpace, upper=upperParametricSapce)
  if(result$convergence != 0)
    result$par = c(NA,NA)
  result$par
}

boostrap <- function(sampledList, estimator, parans){
  sampledListBoot <- replicate(repetitionsBoot, rskellam(n, parans[1], parans[2]), simplify = FALSE)
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
  return (c(real, m, v, b, e))
}

for (par in experiments){
    n <- par[1]
    lambda1 <- par[2]
    lambda2 <- par[3]
    
    cat('Running n=', n, '\n')
      
    start_time <- Sys.time()
    
    set.seed(0)
    
    sampledList <- replicate(repetitions, rskellam(n, lambda1, lambda2), simplify = FALSE)
    
    print("Running mle")
    monteCarlo.mle <- simplify2array(mclapply(sampledList, maximumLikelihoodEstimator, mc.cores = numcpus))
    print("Running mme")
    monteCarlo.mme <- simplify2array(mclapply(sampledList, momentEstimator, mc.cores = numcpus))
    print("Running boot.mme")
    monteCarlo.boot.mme <- simplify2array(mclapply(sampledList, boostrap, momentEstimator, monteCarlo.mme, mc.cores = numcpus))
    print("Running boot.mle")
    monteCarlo.boot.mle <- simplify2array(mclapply(sampledList, boostrap, maximumLikelihoodEstimator, monteCarlo.mle, mc.cores = numcpus))
    
    data <- data.frame(
      mle.y1      = resumeMonteCarlo("y1(mle)     ", monteCarlo.mle[1,], lambda1),
      mle.y2      = resumeMonteCarlo("y2(mle)     ", monteCarlo.mle[2,], lambda2),
      mme.y1      = resumeMonteCarlo("y1(mme)     ", monteCarlo.mme[1,], lambda1),
      mme.y2      = resumeMonteCarlo("y2(mme)     ", monteCarlo.mme[2,], lambda2),
      mle.boot.y1 = resumeMonteCarlo("y1(mle.boot)", monteCarlo.boot.mle[1,], lambda1),
      mle.boot.y2 = resumeMonteCarlo("y2(mle.boot)", monteCarlo.boot.mle[2,], lambda2),
      mme.boot.y1 = resumeMonteCarlo("y1(mme.boot)", monteCarlo.boot.mme[1,], lambda1),
      mme.boot.y2 = resumeMonteCarlo("y2(mme.boot)", monteCarlo.boot.mme[2,], lambda2)
    )
    
    row.names(data) <- c("real", "mean", "variance", "bias", "mse")
    data <- t(data)
    
    name <- paste0("results y1=", lambda1, " y2=",lambda2," n=",n," r=",repetitions," b=", repetitionsBoot ,".csv")
    name <- paste0("./Documents/Skellam_Distribution/Skellam-Distribution/results/", name)
    print(data, digits = 4)
    
    write.csv(data, file = name)
    
    end_time <- Sys.time()
    cat("Elapsed Time:", end_time - start_time, "\n")
  
}
