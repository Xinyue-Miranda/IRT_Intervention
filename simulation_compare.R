library(MASS)
library(mvtnorm)
library(coda)
library(tmvtnorm)
library(nimble)
library(MCMCpack)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(pscl)
library(tidyverse)


##### ---------------------- Simulation ------------------------ #####
#
# Simulating AB Response data for given parameters 
# (theta, lambda and b)
#
sourceCpp("Code/func_0919.cpp")

N    <- 1000
K    <- 100
S    <- 5000
burn <- 2000

mu <- c(0, 0, 0)
I  <- diag(3)
d  <- 3

# parameters
# hierachical
S_0 = I
sigma_sqr <- 1
Sigma     <- I * 100

true_theta  <- rmvnorm(n = N, mean = mu, sigma = 80 * I)
true_lambda <- rmvnorm(n = K, mean = mu, sigma = I)
true_b      <- rnorm(  n = K, mean = 0,  sd    = 1)

# M matrix
M = array(NA, c(3,3,K))
for (k in 1:K) {
  di = rbinom(3, 1, 0.5)
  while(!1 %in% di | !0 %in% di){
    di = rbinom(3, 1, 0.5)
  }
  M[,,k] = diag(di)
}

## Z and Y
true_Z = matrix(NA, N, K)
Y = matrix(NA, N, K)
for (i in 1:N) {
  for (k in 1:K) {
    mu_ik = t(true_lambda[k, ]) %*% M[,,k] %*% true_theta[i, ] - true_b[k]
    true_Z[i, k] = rnorm(1, mean = mu_ik, sd = 1)
    Y[i, k] = ifelse(true_Z[i,k] > 0, 1, 0)
  }
}

# directions
direction <- ifelse(true_lambda > 0, 1, -1)

# Create Fake Points for Y and theta
Flags      <- c(1, -1)
Directions <- c(1, 2, 3)

points <- c()
for (dir in Directions) {
  for (flag in Flags) {
    points <- c(points, ifelse( direction[ , dir] == flag, 1, 0 ))
  }
}

fake_points <- data.frame(
  matrix(points, nrow = 6, ncol = ncol(Y), byrow = TRUE), 
  row.names = paste0("fake_", 1:6)
)
colnames(fake_points) <- colnames(Y)

# Combine with original data
Y <- rbind(fake_points %>% as.matrix(), Y)
cat(c("\nDimension of Y :\n", dim(Y)))

# Create corresponding constraints on thetas
A <- 10
x <- list(
  "1" = c(A,  0,  0),
  "2" = c(-A, 0,  0), 
  "3" = c(0,  A,  0),
  "4" = c(0, -A,  0),
  "5" = c(0,  0,  A),
  "6" = c(0,  0, -A)
)
ind <- 1:6

true_theta <- rbind( data.frame(x) %>% t() %>% as.matrix(), true_theta )

mse <- function(ytrue, ypred){
  ytrue_sd <- apply(ytrue, 2, sd)
  ypred_sd <- apply(ypred, 2, sd)
  ret <- colMeans((ytrue/ytrue_sd-ypred/ypred_sd)^2)
  return(ret)
}

ITER <- 20
MSE <- array(NA, c(ITER, 3, 3))

# Simulation 
for(iter in 1:ITER){
  
  cat(paste0("******************Iter:", iter, "*******************","\n\n"))
  
  
  ##### ------------------------ PCA ----------------------------- #####
  # 
  # Taking PCA on simulated data and
  # keeping the first-3 principal components
  #
  
  p <- prcomp(Y)
  theta_pca <- p$x[,1:3]
  mse1 <- mse(true_theta, theta_pca)
  MSE[iter,,1] <- mse1
  
  ##### ------------------------ pscl ------------------------------ #####
  #
  # Comparing results on simulated data
  #
  
  data_pscl <- rollcall(Y, yea = 1, nay = 0, 
                        legis.names = as.character(1:(N+6)))
  cl <- constrain.legis(data_pscl, x=x, d=d)
  id2Constrained <- ideal(data_pscl,
                          d=d,
                          priors=cl,      ## priors (w constraints)
                          startvals=cl,   ## start value (w constraints)
                          store.item=TRUE,
                          maxiter=5000,
                          burnin=2000,
                          thin=1)
  s <- summary(id2Constrained)
  theta_pscl <- data.frame(s$xm)
  mse2 <- mse(true_theta, theta_pscl)
  MSE[iter,,2] <- mse2
  
  ##### ------------------------ irt ------------------------------ #####
  
  N <- 1006
  
  ## MCMC Setup 
  Z      <- matrix(NA, N, K)
  b      <- rep(0, K)
  theta  <- matrix(0, N, 3)
  lambda <- matrix(0, K, 3)
  
  ## initialize
  inv_v <- array(0.01 * diag(d), c(d, d, N))
  for(i in 1:6){
    # set constraints
    num          <- ind[i]
    theta[num,]  <- x[[i]]
    inv_v[,,num] <- 1e12 * diag(d)
  }
  theta_start <- theta
  
  ## store values
  THETA  <- array(NA, c(N, 3, S-burn))
  LAMBDA <- array(NA, c(K, 3, S-burn))
  
  ## timing
  t_start      <- proc.time()[3]
  t_iter       <- proc.time()[3]
  t_iter_start <- proc.time()[3]
  
  ## updating
  for(s in 1:S){
    ## Timing
    t_total      <- proc.time()[3] - t_start
    t_ETA        <- t_iter * (S - s)
    t_iter       <- proc.time()[3] - t_iter_start
    t_iter_start <- proc.time()[3]
    if (s %% 100 == 0){
      cat(paste("\nIteration:", s, "Last itearation:", round(t_iter, 2), 
                "seconds, Total Elapsed:", round(t_total, 1), 
                "seconds, ETA:", round(t_ETA, 1), "seconds\n")) 
    }
    
    ## update
    Z           <- update_Z(N, K, lambda, M, theta, b, Y)
    b           <- update_b(N, K, sigma_sqr, theta, M, lambda, Z)
    lambda      <- update_lambda(K, theta, M, b, Sigma, Z, d)
    Sigma       <- riwish(3 + K, solve(t(lambda) %*% lambda + S_0))
    sigma_sqr   <- rinvgamma(n = 1, shape = (K + 1)/2, scale = (t(b) %*% b + 1)/2)
    lambda_star <- update_lambda_star(K, lambda, M, d)
    theta       <- update_theta(N, lambda_star, b, Z, inv_v, ind, theta_start, d)
    
    ## storing trace
    if(s > burn) {
      THETA [ , ,s-burn] <- theta
      LAMBDA[ , ,s-burn] <- lambda
    }
  }
  
  theta_irt <- apply(THETA, c(1,2), mean)
  mse3 <- mse(true_theta, theta_irt)
  MSE[iter,,3] <- mse3
  
  cat("\n")
  print(mse1)
  print(mse2)
  print(mse3)
  
  ## Generate new Z and Y
  N = 1000
  for (i in 1:N) {
    for (k in 1:K) {
      mu_ik = t(true_lambda[k, ]) %*% M[,,k] %*% true_theta[i, ] - true_b[k]
      true_Z[i, k] = rnorm(1, mean = mu_ik, sd = 1)
      Y[i, k] = ifelse(true_Z[i,k] > 0, 1, 0)
    }
  }
  Y[1:6,] <- fake_points %>% as.matrix()
  
}

save(MSE,  file = "Results/MSE.RData") 

t(apply(MSE, c(2,3), mean))


