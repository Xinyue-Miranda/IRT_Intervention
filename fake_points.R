## Load Packages
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(tidyverse)
library(coda)
library(MCMCpack)


# ------------------------ MCMC ------------------------------ #

rm(list = ls())

sourceCpp("Code/func_0919.cpp")

country <- "South_Africa"
round   <- "R2" 
S       <- 8000 
burn    <- 3000

# Read Data
Mdata_path <- paste0("Data/M_matrix/", country,"/", country, "_", round, "_M.csv")
M_path     <- paste0("Data/Clean_Data/", country, "/", round, "_M_", country, ".rds")
Y_path     <- paste0("Data/Clean_Data/", country, "/", round, "_Y_", country, ".rds")

Mdata <- read.csv(Mdata_path, fileEncoding = "UTF-8-BOM")
M     <- readRDS(M_path)
Y     <- readRDS(Y_path)


# merge on Question Numbers
Mdata <- Mdata %>%
  dplyr::select(-8) %>%
  rename(
    question  = Survey.Question,
    theta_1   = Theta.1,
    theta_1_d = Theta.1.Direction,
    theta_2   = Theta.2,
    theta_2_d = Theta.2.Direction,
    theta_3   = Theta.3,
    theta_3_d = Theta.3.Direction
  ) %>%
  filter(
    question %in% colnames(Y)
  )

direction <- Mdata[str_detect(colnames(Mdata), "d")]



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
Y <- rbind(fake_points, Y) %>% as.matrix()
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


# MCMC
## Setup
#set.seed(20190717)
N    <- nrow(Y)
K    <- ncol(Y)

d         <- 3
I         <- diag(d)
S_0       <- I
sigma_sqr <- 1
Sigma     <- I * 100


## No Simulation for Real Data
M <- M[,,1:K]
Y <- Y[1:N, 1:K]


## MCMC Setup 
Z      <- matrix(NA, N, K)
b      <- rep(0, K)
theta  <- matrix(0, N, 3)
lambda <- matrix(0, K, 3)


## initialize
inv_v <- array(0.01 * diag(d), c(d, d, N))
for(i in 1:d+1){
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
              "seconds, ETA:", round(t_ETA, 1), "seconds")) 
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


# --------------------- Diagnostic  ----------------------- #

diagnostic_theta <- function(num, N, K, round, country){
  mcmc_obj <- mcmc(as.data.frame(t(THETA[num, , ])))
  print(summary(mcmc_obj))
  print("Effective Sample Size:")
  print(effectiveSize(mcmc_obj))
  png(paste0("Results/", country, "/", country, "_", round,"_THETA_", num, ".png"))
  plot(mcmc_obj)
  dev.off()
}

diagnostic_lambda <- function(num, N, K, round, country){
  mcmc_obj <- mcmc(as.data.frame(t(LAMBDA[num, , ])))
  print(summary(mcmc_obj))
  print("Effective Sample Size:")
  print(effectiveSize(mcmc_obj))
  png(paste0("Results/", country, "/", country, "_", round,"_LAMBDA_", num, ".png"))
  plot(mcmc_obj)
  dev.off()
}


rand_n <- sample(1:N, 2, replace = FALSE)
rand_k <- sample(1:K, 2, replace = FALSE)
print( paste("Random Ns :", rand_n[1],  rand_n[2]) )
print( paste("Random Ks :", rand_k[1],  rand_k[2]) )

diagnostic_theta (rand_n[1], N, K, round, country)
diagnostic_theta (rand_n[2], N, K, round, country)
diagnostic_lambda(rand_k[1], N, K, round, country)
diagnostic_lambda(rand_k[2], N, K, round, country)


## save samplers
save(THETA,  file = paste0("Results/",country, "/", country, "_", round, "_THETA.RData")) 
save(LAMBDA, file = paste0("Results/",country, "/", country, "_", round, "_LAMBDA.RData")) 


