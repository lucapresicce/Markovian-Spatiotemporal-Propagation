# -------------------------------------------------------------------
# AMORTIZED BAYESIAN FORECAST - GENERATIVE FUNCTIONS
# -------------------------------------------------------------------

# Working packages
library(mniw)
library(classInt)
library(RColorBrewer)
library(sp)
library(fields)
library(parallel)
library(foreach)
library(doParallel)
library(akima)
library(ggplot2)

# Associated R package
# devtools::install_github("lucapresicce/spFFBS") # if not previously installed
library(spFFBS)
library(spBPS)

# Data generation ---------------------------------------------------------

tendency_gen <- function(n = 50, q = 2, p = 2, tmax = 10, phi, seed = 42) {
  
  set.seed(42)
  # p <- 2
  # q <- 2
  fixed_crd <- matrix(runif(n * 2), ncol = 2)
  Sigma <- matrix(c(1, -0.3, -0.3, 1), q, q)
  Rphi_val <- exp(-phi * arma_dist(fixed_crd))
  fixed_W <- mniw::rMNorm(1, Lambda = matrix(0, n, q), SigmaR = Rphi_val, SigmaC = Sigma)
  
  
  set.seed(seed)
  fixed_crd <- matrix(runif(n * 2), ncol = 2)
  
  D   <- spBPS::arma_dist(fixed_crd)
  K   <- exp(-phi*D)
  W <- rbind( cbind(diag(p), matrix(0, p, n)), cbind(matrix(0, n, p), K) )
  
  # Prior information
  m0     <- matrix(0, n+p, q)
  C0 <- rbind( cbind(diag(0.005, p), matrix(0, p, n)), cbind(matrix(0, n, p), K) )
  nu0 <- p+2
  Psi0 <- diag(q)
  
  # Initial state
  theta0 <- mniw::rMNorm(n = 1, Lambda = m0, SigmaR = C0, SigmaC = Sigma)
  
  # Generate dynamic parameters from the model
  G <- diag(n+p)
  theta <- array(0, c(n+p, q, tmax+1))
  set.seed(seed)
  for (t in 1:(tmax+1)) {
    if (t >= 2) {  
      theta[,,t] <- G %*% theta[,,t-1] + mniw::rMNorm(n = 1, Lambda = m0, SigmaR = W, SigmaC = Sigma)
    } else {
      theta[,,t] <- G %*% theta0 + mniw::rMNorm(n = 1, Lambda = m0, SigmaR = W, SigmaC = Sigma)
    }
  }
  
  return(list("Theta" = theta, "crd" = fixed_crd, "G" = G))
  
}

data_gen <- function(n = 50, q = 2, p = 2, tmax = 10, Theta)  {
  
  # p <- 2
  # q <- 2
  Sigma <- matrix(c(1, -0.3, -0.3, 1), q, q)
  tau <- 0.8
  a <- ((1/tau)-1)
  V <- a*diag(n)
  
  X     <- cbind(rep(1, n), matrix(runif((p - 1) * n, -1, 1), ncol = p - 1))
  P     <- cbind(X, diag(n))
  Y <- array(0, c(n, q, tmax+1))
  for (t in 1:(tmax+1)) {
    if (t >= 2) {  
      
      Y[,,t]     <- P %*% Theta[,,t] + mniw::rMNorm(n = 1, Lambda = matrix(0, n, q), SigmaR = V, SigmaC = Sigma)
    } else {
      
      Y[,,t]     <- P %*% Theta[,,t] + mniw::rMNorm(n = 1, Lambda = matrix(0, n, q), SigmaR = V, SigmaC = Sigma)
    }
  }
  
  X     <- replicate(expr = X, n = tmax+1)
  P     <- abind::abind(X, replicate(expr = diag(n), n = tmax+1), along = 2)
  
  return(abind::abind(Y, P, along=2))
}

# Calculate post ----------------------------------------------------------

calculate_post <- function(Z, fixed, q = 2, p = 2) {
  
  # dimensions
  # p <- 2
  # q <- 2
  n <- dim(Z)[1]
  tmax <- dim(Z)[3]-1
  
  #  Sampling data
  Y    <- Z[, 1:q, 1:tmax]
  X    <- Z[, (q+1):(q+p),]
  P    <- Z[,-(1:q),]
  G     <- fixed$G
  crd   <- fixed$crd
  D   <- spBPS::arma_dist(crd)
  K   <- exp(-D)
  
  # Prior information
  m0     <- matrix(0, n+p, q)
  C0 <- rbind( cbind(diag(0.005, p), matrix(0, p, n)), cbind(matrix(0, n, p), K) )
  nu0 <- p+2
  Psi0 <- diag(q)
  prior <- list("m" = m0, "C" = C0, "nu" = nu0, "Psi" = Psi0)
  
  # setting up DYNBPS
  tau_seq <- c(0.7, 0.8, 0.9)
  phi_seq <- c(2, 4, 6)
  par_grid <- spBPS::expand_grid_cpp(rev(tau_seq), rev(phi_seq))
  J <- nrow(par_grid)
  
  # results
  G3 <- array(G, c(dim(G)[1], dim(G)[2], dim(P)[3]))
  out <- spFFBS::spFFBS(Y = Y, G = G3, P = P, D = D,
                        grid = list(tau = tau_seq, phi = phi_seq), prior = prior, L = 200, tnew = 1, do_forecast = T,
                        do_spatial = T, spatial = list(crd = crd, crdtilde = crd, Xtilde = X, t = tmax))
  
  L <- 200
  Ysp_pred   <- sapply(1:L, function(l){out$spatial[[1]][[l]][1:n,]}   , simplify = "array")
  Omega_pred <- sapply(1:L, function(l){out$spatial[[1]][[l]][-(1:n),]}, simplify = "array")
  
  # RMSPE
  Ysp_postmean <- apply(Ysp_pred, c(1,2), mean)
  Ysp_upp      <- apply(Ysp_pred, c(1,2), quantile, 0.975)
  Ysp_low      <- apply(Ysp_pred, c(1,2), quantile, 0.025)
  # colMeans(sqrt((Z[,1:q,tmax+1] - Ysp_postmean)^2))
  
  Omega_postmean <- apply(Omega_pred, c(1,2), mean)
  Omega_upp      <- apply(Omega_pred, c(1,2), quantile, 0.975)
  Omega_low      <- apply(Omega_pred, c(1,2), quantile, 0.025)
  # colMeans(sqrt((fixed$Theta[-(1:p),,tmax+1] - Omega_postmean)^2))
  
  # return posterior quantiles
  Ups_50  <- cbind(Ysp_postmean, Omega_postmean)
  Ups_025 <- cbind(Ysp_low, Omega_low)
  Ups_975 <- cbind(Ysp_upp, Omega_upp)
  W <- abind::abind(Ups_50, Ups_025, Ups_975, along = 3)
  
  return(W)
  
}


# Plot --------------------------------------------------------------------

# Plot surface function
plot_surface_interp <- function(mat, fixed_crd, title = NULL, component = 1, grid_res = 100) {
  library(RColorBrewer)
  vals <- mat[, component]
  interp_result <- with(data.frame(x = fixed_crd[, 1], y = fixed_crd[, 2], z = vals),
                        interp(x, y, z, nx = grid_res, ny = grid_res, duplicate = "mean"))
  df_grid <- data.frame(
    x = rep(interp_result$x, times = length(interp_result$y)),
    y = rep(interp_result$y, each = length(interp_result$x)),
    value = as.vector(interp_result$z)
  )
  
  ggplot(df_grid, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdBu"))) + 
    coord_equal() +
    ggtitle(title) +
    theme_minimal()
}

# Check -------------------------------------------------------------------

# compute
fixed <- tendency_gen(n = 100, q = 2, p = 2, tmax = 10, phi = 4, seed = 1997)
Z     <- data_gen(n = 100, tmax = 10, Theta = fixed$Theta)
W     <- calculate_post(Z = Z, fixed = fixed)

# Y
cowplot::plot_grid(
  plot_surface_interp(mat = Z[,1:q,tmax+1], fixed_crd = fixed$crd, component = 1, title = "TRUE Y1"),
  plot_surface_interp(mat = Z[,1:q,tmax+1], fixed_crd = fixed$crd, component = 2, title = "TRUE Y2"),
  plot_surface_interp(mat = W[,1:q,1], fixed_crd = fixed$crd, component = 1     , title = "PRED Y1"),
  plot_surface_interp(mat = W[,1:q,1], fixed_crd = fixed$crd, component = 2     , title = "PRED Y2"),
  nrow = 2)

# Omega
cowplot::plot_grid(
  plot_surface_interp(mat = fixed$Theta[-(1:p),,tmax+1], fixed_crd = fixed$crd, component = 1, title = "TRUE OM1"),
  plot_surface_interp(mat = fixed$Theta[-(1:p),,tmax+1], fixed_crd = fixed$crd, component = 2, title = "TRUE OM2"),
  plot_surface_interp(mat = W[,-(1:q),1], fixed_crd = fixed$crd, component = 1               , title = "PRED OM1"),
  plot_surface_interp(mat = W[,-(1:q),1], fixed_crd = fixed$crd, component = 2               , title = "PRED OM2"),
  nrow = 2)
