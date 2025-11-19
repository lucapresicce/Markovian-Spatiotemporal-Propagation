### M-CLOSED/M-OPEN SETTINGS - SIMULATION #################################

rm(list = ls())
gc()

# Packages ----------------------------------------------------------------

# Working packages
library(mniw)
library(classInt)
library(RColorBrewer)
library(sp)
library(fields)
library(parallel)
library(foreach)
library(doParallel)

# Associated R package
# devtools::install_github("lucapresicce/spFFBS") # if not previously installed
library(spFFBS)
library(spBPS)

# Replications ------------------------------------------------------------

# Synthetic data dimensions
tmax <- 20
tnew <- 5
n    <- 100
u    <- 50
q    <- 3
p    <- 2

# parameters
Sigma  <- matrix(c(1, -0.3, 0.6, -0.3, 1.2, 0.4, 0.6, 0.4, 1), q, q)
phi <- 8
tau <- 0.80
a <- ((1/tau)-1)
V <- a*diag((n+u))

set.seed(97)
# generate constant data structure
coords <- matrix(runif((n+u) * 2), ncol = 2)
D <- arma_dist(coords); gc()
K <- exp(-phi * D)
W <- rbind( cbind(diag(p), matrix(0, p, (n+u))), cbind(matrix(0, (n+u), p), K) )
rm("D"); gc()

# Prior information
m0     <- matrix(0, (n+u)+p, q)
C0 <- rbind( cbind(diag(0.005, p), matrix(0, p, (n+u))), cbind(matrix(0, (n+u), p), K) )

# Initial state
theta0 <- mniw::rMNorm(n = 1, Lambda = m0, SigmaR = C0, SigmaC = Sigma)

# Generate dynamic parameters from the model
G <- array(0, c((n+u)+p, (n+u)+p, tmax+tnew))
theta <- array(0, c((n+u)+p, q, tmax+tnew))
set.seed(4-8-15-16-23-42)
for (t in 1:(tmax+tnew)) {
  if (t >= 2) {  
    G[,,t]     <- diag(p+n+u)
    theta[,,t] <- G[,,t] %*% theta[,,t-1] + mniw::rMNorm(n = 1, Lambda = m0, SigmaR = W, SigmaC = Sigma)
  } else {
    G[,,t]     <- diag(p+n+u)
    theta[,,t] <- G[,,t] %*% theta0 + mniw::rMNorm(n = 1, Lambda = m0, SigmaR = W, SigmaC = Sigma)
  }
}


# replications loop
BB <- 5
out_list <- vector(mode = "list", length = BB)
pb <- txtProgressBar(min = 0, max = BB, style = 3)
for (b in 1:BB) {
  
  cat("Computation started ... \n")
  
  # define models sub-list
  replication <- vector(mode = "list", length = 4)
  
  # Fix the seed for reproducibility
  set.seed(b)
  
  # Data generation -------------------------------------------------------
  
  # fixed-time covariates
  # Generate data from the model
  V <- a*diag((n+u))
  X <- array(0, c((n+u), p, tmax+tnew))
  P <- array(0, c((n+u), (n+u)+p, tmax+tnew))
  Y <- array(0, c((n+u), q, tmax+tnew))
  # set.seed(4-8-15-16-23-42)
  for (t in 1:(tmax+tnew)) {
    if (t >= 2) {  
      
      X[,,t]     <- matrix(runif((n+u)*p), (n+u), p)
      P[,,t]     <- cbind(X[,,t], diag((n+u)))
      Y[,,t]     <- P[,,t] %*% theta[,,t] + mniw::rMNorm(n = 1, Lambda = matrix(0, (n+u), q), SigmaR = V, SigmaC = Sigma)
    } else {
      
      X[,,t]     <- matrix(runif((n+u)*p), (n+u), p)
      P[,,t]     <- cbind(X[,,t], diag((n+u)))
      Y[,,t]     <- P[,,t] %*% theta[,,t] + mniw::rMNorm(n = 1, Lambda = matrix(0, (n+u), q), SigmaR = V, SigmaC = Sigma)
    }
  }
  
  # save test data
  replication[[4]] <- Y
  
  # Unobserved locations data
  Ytilde <- Y[-(1:n),,]
  thetatilde <- theta[-(1:(n+p)),,]
  Xtilde <- X[-(1:n),,]
  crdtilde <- coords[-(1:n),]
  Dtilde   <- as.matrix(dist(crdtilde))
  Ktilde   <- exp(-phi*Dtilde)
  
  #  Sampling data
  Y <- Y[(1:n),,1:tmax]
  X <- X[(1:n),,]
  P     <- P[(1:n),1:(n+p),]
  G     <- G[1:(n+p), 1:(n+p),]
  crd <- coords[1:n,]
  D   <- as.matrix(dist(crd))
  K   <- exp(-D)
  W <- rbind( cbind(diag(p), matrix(0, p, n)), cbind(matrix(0, n, p), K) )
  V <- a*diag((n))
  
  # Prior information
  m0     <- matrix(0, n+p, q)
  C0 <- rbind( cbind(diag(0.005, p), matrix(0, p, n)), cbind(matrix(0, n, p), K) )
  nu0 <- 3
  Psi0 <- diag(q)
  
  # prior information
  prior <- list("m" = m0, "C" = C0, "nu" = nu0, "Psi" = Psi0)
  
  # free memory
  gc()
  
  cat("Data generated ... \n")
  
  # True Model - Sequential Learning --------------------------------------
  
  # Define competetive (J) models
  tau_seq <- tau
  phi_seq <- phi
  par_grid <- spBPS::expand_grid_cpp(rev(tau_seq), rev(phi_seq))
  J <- nrow(par_grid)
  
  # posterior samples and random time points for interpolations
  L <- 200
  set.seed(b)
  t <- sample(1:(tmax-1), 1, F)
  to <- sample((tmax+1):(tmax+tnew), 1, F)

  out <- spFFBS::spFFBS(Y = Y, G = G, P = P, D = D,
                        grid = list(tau = tau_seq, phi = phi_seq), prior = prior, L = L,
                        do_BS = T, do_forecast = T, do_spatial = T, tnew = tnew,
                        spatial = list(crd = crd, crdtilde = crdtilde, Xtilde = Xtilde,
                                       t = c(t, to)))
  # Global weights
  Wglobal <- out$Wglobal
  
  # collect Sigma posterior inference
  indL <- sample(1:J, L, Wglobal, rep = T)
  Sigma_post <- sapply(1:L, function(l) {
    mniw::riwish(1, nu = out$FF[[tmax]]$filtered_results[[indL[l]]]$nu,
                 Psi = out$FF[[tmax]]$filtered_results[[indL[l]]]$Psi) },
    simplify = "array")
  
  # collect Theta posterior inference
  theta_post <- sapply(1:tmax, function(t){ out$BS[[t]] }, simplify = "array")
  
  # temporal forecasts
  Y_pred <- out$forecast$Y_pred
  
  # spatial interpolations
  out_SI_i <- out$spatial[[1]]
  out_SI_o <- out$spatial[[2]]
  
  # Replicazione
  replication[[1]] <- list("Sigma posterior" = Sigma_post,
                           "theta posterior" = theta_post,
                           "prediction"   = Y_pred,
                           "interpol in"  = out_SI_i,
                           "interpol out" = out_SI_o,
                           "weights" = Wglobal)
  
  cat("True Model completed ... \n")
  
  # DYNBPS closed - Parallel Learning --------------------------------------
  
  # Define competetive (J) models
  tau_seq <- c(0.7, 0.8, 0.9)
  phi_seq <- c(6, 8, 10)
  par_grid <- spBPS::expand_grid_cpp(rev(tau_seq), rev(phi_seq))
  J <- nrow(par_grid)
  
  # posterior samples and random time points for interpolations
  L <- 200
  set.seed(b)
  t <- sample(1:(tmax-1), 1, F)
  to <- sample((tmax+1):(tmax+tnew), 1, F)
  
  out <- spFFBS::spFFBS(Y = Y, G = G, P = P, D = D,
                        grid = list(tau = tau_seq, phi = phi_seq), prior = prior, L = L,
                        do_BS = T, do_forecast = T, do_spatial = T, tnew = tnew,
                        spatial = list(crd = crd, crdtilde = crdtilde, Xtilde = Xtilde,
                                       t = c(t, to)))
  # Global weights
  Wglobal <- out$Wglobal
  
  # collect Sigma posterior inference
  indL <- sample(1:J, L, Wglobal, rep = T)
  Sigma_post <- sapply(1:L, function(l) {
    mniw::riwish(1, nu = out$FF[[tmax]]$filtered_results[[indL[l]]]$nu,
                 Psi = out$FF[[tmax]]$filtered_results[[indL[l]]]$Psi) },
    simplify = "array")
  
  # collect Theta posterior inference
  theta_post <- sapply(1:tmax, function(t){ out$BS[[t]] }, simplify = "array")
  
  # temporal forecasts
  Y_pred <- out$forecast$Y_pred
  
  # spatial interpolations
  out_SI_i <- out$spatial[[1]]
  out_SI_o <- out$spatial[[2]]
  
  # Replicazione
  replication[[2]] <- list("Sigma posterior" = Sigma_post,
                           "theta posterior" = theta_post,
                           "prediction"   = Y_pred,
                           "interpol in"  = out_SI_i,
                           "interpol out" = out_SI_o,
                           "weights" = Wglobal)
  
  cat("DYNBPS completed - closed ... \n")
  
  # DYNBPS open - Parallel Learning ----------------------------------------
  
  # Define competetive (J) models
  # set.seed(97)
  (tau_seq <- runif(3, 0.5, 1))
  # set.seed(97)
  (phi_seq <- sample(1:50, 3, F))
  par_grid <- spBPS::expand_grid_cpp(rev(tau_seq), rev(phi_seq))
  J <- nrow(par_grid)
  
  # posterior samples and random time points for interpolations
  L <- 200
  set.seed(b)
  t <- sample(1:(tmax-1), 1, F)
  to <- sample((tmax+1):(tmax+tnew), 1, F)
  
  out <- spFFBS::spFFBS(Y = Y, G = G, P = P, D = D,
                        grid = list(tau = tau_seq, phi = phi_seq), prior = prior, L = L,
                        do_BS = T, do_forecast = T, do_spatial = T, tnew = tnew,
                        spatial = list(crd = crd, crdtilde = crdtilde, Xtilde = Xtilde,
                                       t = c(t, to)))
  # Global weights
  Wglobal <- out$Wglobal
  
  # collect Sigma posterior inference
  indL <- sample(1:J, L, Wglobal, rep = T)
  Sigma_post <- sapply(1:L, function(l) {
    mniw::riwish(1, nu = out$FF[[tmax]]$filtered_results[[indL[l]]]$nu,
                 Psi = out$FF[[tmax]]$filtered_results[[indL[l]]]$Psi) },
    simplify = "array")
  
  # collect Theta posterior inference
  theta_post <- sapply(1:tmax, function(t){ out$BS[[t]] }, simplify = "array")
  
  # temporal forecasts
  Y_pred <- out$forecast$Y_pred
  
  # spatial interpolations
  out_SI_i <- out$spatial[[1]]
  out_SI_o <- out$spatial[[2]]
  
  # Replicazione
  replication[[3]] <- list("Sigma posterior" = Sigma_post,
                           "theta posterior" = theta_post,
                           "prediction"   = Y_pred,
                           "interpol in"  = out_SI_i,
                           "interpol out" = out_SI_o,
                           "weights" = Wglobal)
  
  cat("DYNBPS completed - open ... \n")
  
  # assign models sub-list ------------------------------------------------
  out_list[[b]] <- replication
  
  # free memory
  rm(list = c("replication"))
  
  # pb
  setTxtProgressBar(pb, b)
  
  
}

# free the memory
gc()

# save the replications
# save(out_list, file = "replications_results.RData")
# load("replications_results.RData")

