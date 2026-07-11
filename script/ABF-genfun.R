#####################################################################
# AMORTIZED BAYESIAN FORECAST - GENERATIVE FUNCTIONS
#####################################################################

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
library(abind)
library(MBA)

# associated packages
library(spFFBS)
library(spBPS)

# new C++ routines and wrapper
Rcpp::sourceCpp("script/FFBS-DYNBPS-struct-v2.cpp")
source("script/FFBS-DYNBPS-v2.R")


# Data generation ---------------------------------------------------------

tendency_gen <- function(n = 50, q = 2, p = 2, tmax = 10, phi, seed = 42) {

  set.seed(42)
  fixed_crd <- matrix(runif(n * 2), ncol = 2)
  Sigma     <- matrix(c(1, -0.3, -0.3, 1), q, q)
  Rphi_val  <- exp(-phi * spBPS:::arma_dist(fixed_crd))
  fixed_W   <- mniw::rMNorm(1, Lambda = matrix(0, n, q),
                             SigmaR = Rphi_val, SigmaC = Sigma)

  set.seed(seed)
  fixed_crd <- matrix(runif(n * 2), ncol = 2)

  D  <- spBPS:::arma_dist(fixed_crd)
  K  <- exp(-phi * D)
  W  <- rbind(cbind(diag(p),       matrix(0, p, n)),
              cbind(matrix(0, n, p), K))

  m0  <- matrix(0, n + p, q)
  C0  <- rbind(cbind(diag(0.005, p), matrix(0, p, n)),
               cbind(matrix(0, n, p), K))
  nu0  <- p + 2
  Psi0 <- diag(q)

  theta0 <- mniw::rMNorm(n = 1, Lambda = m0, SigmaR = C0, SigmaC = Sigma)

  G     <- diag(n + p)
  theta <- array(0, c(n + p, q, tmax + 1))
  set.seed(seed)
  for (t in 1:(tmax + 1)) {
    if (t >= 2) {
      theta[,, t] <- G %*% theta[,, t-1] +
        mniw::rMNorm(n = 1, Lambda = m0, SigmaR = W, SigmaC = Sigma)
    } else {
      theta[,, t] <- G %*% theta0 +
        mniw::rMNorm(n = 1, Lambda = m0, SigmaR = W, SigmaC = Sigma)
    }
  }

  return(list("Theta" = theta, "crd" = fixed_crd, "G" = G))
}


data_gen <- function(n = 50, q = 2, p = 2, tmax = 10, Theta) {

  Sigma <- matrix(c(1, -0.3, -0.3, 1), q, q)
  tau   <- 0.8
  a     <- (1 / tau) - 1
  V     <- a * diag(n)

  X <- cbind(rep(1, n), matrix(runif((p - 1) * n, -1, 1), ncol = p - 1))
  P <- cbind(X, diag(n))

  Y <- array(0, c(n, q, tmax + 1))
  for (t in 1:(tmax + 1)) {
    Y[,, t] <- P %*% Theta[,, t] +
      mniw::rMNorm(n = 1, Lambda = matrix(0, n, q), SigmaR = V, SigmaC = Sigma)
  }

  X <- replicate(expr = X, n = tmax + 1)
  P <- abind::abind(X, replicate(expr = diag(n), n = tmax + 1), along = 2)

  return(abind::abind(Y, P, along = 2))
}


# Calculate posterior - DYNBPS teacher ------------------------------------

calculate_post <- function(Z, fixed, q = 2, p = 2, n_threads = 1L) {

  n    <- dim(Z)[1]
  tmax <- dim(Z)[3] - 1L

  Y       <- Z[, 1:q,         1:tmax,        drop = FALSE]  # n x q x tmax
  X_train <- Z[, (q+1):(q+p), 1:tmax,        drop = FALSE]  # n x p x tmax
  X_all   <- Z[, (q+1):(q+p), ,              drop = FALSE]  # n x p x (tmax+1)
  Xtilde  <- Z[, (q+1):(q+p), tmax + 1L,    drop = FALSE]  # n x p x 1

  crd <- fixed$crd
  D   <- spBPS:::arma_dist(crd)

  m0   <- matrix(0, n + p, q)
  C0   <- rbind(cbind(diag(0.005, p), matrix(0, p, n)),
                cbind(matrix(0, n, p), exp(-D)))          # phi=1 for C0
  nu0  <- p + 2L
  Psi0 <- diag(q)
  prior_v2 <- make_prior_v2(m0, C0, nu0 = nu0, Psi0 = Psi0, p = p)

  tau_seq <- c(0.7, 0.8, 0.9)
  phi_seq <- c(2, 4, 6)

  out <- spFFBS_v2(
    Y           = Y,
    X           = X_train,
    D           = D,
    grid        = list(tau = tau_seq, phi = phi_seq),
    prior       = prior_v2,
    G_beta      = NULL,
    rho         = 1.0,
    do_BS       = FALSE,
    do_forecast = FALSE,
    do_spatial  = TRUE,
    L           = 200L,
    tnew        = 1L,
    X_all       = X_all,
    spatial     = list(
      crd      = crd,
      crdtilde = crd,
      Xtilde   = Xtilde,
      t        = as.integer(tmax + 1L)
    ),
    n_threads = as.integer(n_threads),
    verbose   = FALSE
  )

  Ysp_pred   <- out$spatial[[1L]]$Y      # n x q x L
  Omega_pred <- out$spatial[[1L]]$Omega  # n x q x L

  Ups_50  <- cbind(apply(Ysp_pred,   c(1,2), mean),
                   apply(Omega_pred, c(1,2), mean))
  Ups_025 <- cbind(apply(Ysp_pred,   c(1,2), quantile, 0.025),
                   apply(Omega_pred, c(1,2), quantile, 0.025))
  Ups_975 <- cbind(apply(Ysp_pred,   c(1,2), quantile, 0.975),
                   apply(Omega_pred, c(1,2), quantile, 0.975))

  W <- abind::abind(Ups_50, Ups_025, Ups_975, along = 3)
  return(W)
}


# Calculate posterior - MCMC teacher -------------------------------------

calculate_post_mcmc <- function(Z, fixed, q = 2, p = 2,
                                n_iter          = 2000L,
                                n_burnin        = 1000L,
                                L               = 200L,
                                phi_prior_shape = 2,
                                phi_prior_rate  = 0.5) {

  n    <- dim(Z)[1]
  tmax <- dim(Z)[3] - 1L

  Y_train <- Z[, 1:q,         1:tmax, drop = FALSE]   # n x q x tmax
  X_fixed <- Z[, (q+1):(q+p), 1L]                     # n x p (time-invariant)
  P_mat   <- cbind(X_fixed, diag(n))                  # n x (p+n)

  crd   <- fixed$crd
  D_loc <- spBPS:::arma_dist(crd)

  m0   <- matrix(0, n + p, q)
  nu0  <- p + 2L
  Psi0 <- diag(q)

  lgamma_mq_loc <- function(a, q_) sum(lgamma(a - (seq_len(q_) - 1L) / 2))

  kf_logml_diag <- function(tau, phi) {
    if (tau <= 0 || tau >= 1 || phi <= 0) return(-1e10)
    tryCatch({
      a_v    <- (1 - tau) / tau
      V_kf   <- a_v * diag(n)
      W_phi  <- rbind(cbind(diag(p),        matrix(0, p, n)),
                      cbind(matrix(0, n, p), exp(-phi * D_loc)))
      C_kf   <- rbind(cbind(diag(0.005, p), matrix(0, p, n)),
                      cbind(matrix(0, n, p), exp(-phi * D_loc)))
      m_kf   <- m0
      nu_kf  <- max(nu0, q + 2L)
      Psi_kf <- Psi0
      log_ml <- 0
      for (tt in seq_len(tmax)) {
        R_t  <- C_kf + W_phi
        Qt   <- P_mat %*% R_t %*% t(P_mat) + V_kf
        Qc   <- tryCatch(chol(Qt + diag(n) * 1e-10),
                         error = function(e) diag(n))
        Qi   <- chol2inv(Qc)
        et   <- Y_train[,, tt] - P_mat %*% m_kf
        Psi_n <- Psi_kf + t(et) %*% Qi %*% et
        nu_n  <- nu_kf + n
        lQ   <- 2 * sum(log(diag(Qc)))
        lP0  <- as.numeric(determinant(Psi_kf, log = TRUE)$modulus)
        lPn  <- as.numeric(determinant(Psi_n,  log = TRUE)$modulus)
        log_ml <- log_ml + (
          -0.5 * n * q * log(pi) +
            lgamma_mq_loc(nu_n / 2, q) - lgamma_mq_loc(nu_kf / 2, q) +
            0.5 * nu_kf * lP0 - 0.5 * nu_n * lPn - 0.5 * q * lQ)
        Kt     <- R_t %*% t(P_mat) %*% Qi
        m_kf   <- m_kf + Kt %*% et
        C_kf   <- (diag(n + p) - Kt %*% P_mat) %*% R_t
        C_kf   <- (C_kf + t(C_kf)) / 2
        nu_kf  <- nu_n; Psi_kf <- Psi_n
      }
      log_ml
    }, error = function(e) -1e10)
  }

  log_post_mc <- function(lt, lp) {
    kf_logml_diag(plogis(lt), exp(lp)) +
      dgamma(exp(lp), shape = phi_prior_shape, rate = phi_prior_rate, log = TRUE)
  }

  lt_c     <- qlogis(0.75); lp_c <- log(3)
  lp_c_val <- log_post_mc(lt_c, lp_c)
  chain    <- matrix(NA_real_, n_iter, 2L,
                     dimnames = list(NULL, c("tau", "phi")))
  acc <- 0L; sc <- c(0.15, 0.20)

  t0_chain <- proc.time()
  for (it in seq_len(n_iter)) {
    lt_p     <- lt_c + rnorm(1, 0, sc[1])
    lp_p     <- lp_c + rnorm(1, 0, sc[2])
    lp_p_val <- log_post_mc(lt_p, lp_p)
    if (log(runif(1)) < lp_p_val - lp_c_val) {
      lt_c <- lt_p; lp_c <- lp_p; lp_c_val <- lp_p_val; acc <- acc + 1L
    }
    chain[it,] <- c(plogis(lt_c), exp(lp_c))
    if (it %% 200L == 0L && it <= n_burnin) {
      rate <- acc / it
      if (rate > 0.44) sc <- sc * 1.25 else if (rate < 0.20) sc <- sc / 1.25
    }
  }
  time_chain <- (proc.time() - t0_chain)["elapsed"]

  post   <- chain[(n_burnin + 1L):n_iter,, drop = FALSE]
  keep   <- sample(nrow(post), L, replace = FALSE)
  tau_dr <- post[keep, "tau"]
  phi_dr <- post[keep, "phi"]

  Y_pred_mc  <- array(NA_real_, c(n, q, L))
  Om_pred_mc <- array(NA_real_, c(n, q, L))

  t0_pred <- proc.time()
  for (l in seq_len(L)) {
    tau_l <- tau_dr[l]; phi_l <- phi_dr[l]
    a_v_l <- (1 - tau_l) / tau_l
    V_l   <- a_v_l * diag(n)

    W_phi_l <- rbind(cbind(diag(p),        matrix(0, p, n)),
                     cbind(matrix(0, n, p), exp(-phi_l * D_loc)))
    C_kf_l  <- rbind(cbind(diag(0.005, p), matrix(0, p, n)),
                     cbind(matrix(0, n, p), exp(-phi_l * D_loc)))
    m_kf_l  <- m0; nu_kf_l <- max(nu0, q + 2L); Psi_kf_l <- Psi0

    for (tt in seq_len(tmax)) {
      R_t_l  <- C_kf_l + W_phi_l
      Qt     <- P_mat %*% R_t_l %*% t(P_mat) + V_l
      Qc     <- tryCatch(chol(Qt + diag(n) * 1e-10), error = function(e) diag(n))
      Qi     <- chol2inv(Qc)
      et     <- Y_train[,, tt] - P_mat %*% m_kf_l
      Psi_kf_l <- Psi_kf_l + t(et) %*% Qi %*% et
      nu_kf_l  <- nu_kf_l + n
      Kt     <- R_t_l %*% t(P_mat) %*% Qi
      m_kf_l <- m_kf_l + Kt %*% et
      C_kf_l <- (diag(n + p) - Kt %*% P_mat) %*% R_t_l
      C_kf_l <- (C_kf_l + t(C_kf_l)) / 2
    }

    Sig_l <- tryCatch(mniw::riwish(1L, Psi = Psi_kf_l, nu = nu_kf_l),
                      error = function(e) diag(q))
    LcS   <- tryCatch(t(chol(Sig_l + diag(q) * 1e-8)),
                      error = function(e) diag(q))

    # One-step-ahead state prior: R_{T+1} = C_T + W_phi
    R_tp1 <- C_kf_l + W_phi_l
    R_tp1 <- (R_tp1 + t(R_tp1)) / 2
    LcR   <- tryCatch(t(chol(R_tp1 + diag(n + p) * 1e-8)),
                      error = function(e) diag(n + p))

    Theta_tp1 <- m_kf_l +
      LcR %*% matrix(rnorm((n + p) * q), n + p, q) %*% LcS

    Y_tp1_l <- P_mat %*% Theta_tp1 +
      sqrt(a_v_l) * diag(n) %*% matrix(rnorm(n * q), n, q) %*% LcS

    Y_pred_mc [,, l] <- Y_tp1_l
    Om_pred_mc[,, l] <- Theta_tp1[(p + 1L):(p + n),, drop = FALSE]
  }
  time_pred <- (proc.time() - t0_pred)["elapsed"]

  Ups_50  <- cbind(apply(Y_pred_mc,  c(1,2), mean),
                   apply(Om_pred_mc, c(1,2), mean))
  Ups_025 <- cbind(apply(Y_pred_mc,  c(1,2), quantile, 0.025),
                   apply(Om_pred_mc, c(1,2), quantile, 0.025))
  Ups_975 <- cbind(apply(Y_pred_mc,  c(1,2), quantile, 0.975),
                   apply(Om_pred_mc, c(1,2), quantile, 0.975))

  W <- abind::abind(Ups_50, Ups_025, Ups_975, along = 3)

  attr(W, "time_chain") <- time_chain
  attr(W, "time_pred")  <- time_pred
  attr(W, "acc_rate")   <- acc / n_iter

  return(W)
}


# Plot surface ------------------------------------------------------------

theme_paper <- function(base = 14, leg = "none") {
  theme_bw(base_size = base) %+replace% theme(
    title             = element_text(size = base + 2, face = "bold"), 
    strip.text        = element_text(size = base, face = "bold"),
    strip.background  = element_rect(fill = "#eceff1", colour = "grey70"),
    # axis.text         = element_text(size = base - 2, face = "bold"),
    axis.text.x       = element_text(size = base - 2, hjust = 0.5),
    axis.text.y       = element_text(size = base - 2, hjust = 0.5),
    axis.title        = element_text(size = base, face = "bold"),
    axis.title.x      = element_text(size = base, face = "bold"),
    axis.title.y      = element_text(size = base, face = "bold", angle = 90),
    legend.text       = element_text(size = base + 2, face = "bold"),
    legend.title      = element_blank(),
    legend.key.size   = unit(1.3, "lines"),
    legend.position   = leg,
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_blank()
  )
}

plot_surface_interp <- function(mat, fixed_crd, title = NULL,
                                component = 1,
                                grid_res = 200,
                                limits = NULL) {

  vals <- mat[, component]

  x_seq <- seq(0, 1, length.out = grid_res)
  y_seq <- seq(0, 1, length.out = grid_res)
  
  interp_result <- with(
    data.frame(x = fixed_crd[, 1],
               y = fixed_crd[, 2],
               z = vals),
    interp(x, y, z,
           xo = x_seq,
           yo = y_seq,
           duplicate = "mean")
  )

  df_grid <- data.frame(
    x     = rep(interp_result$x, times = length(interp_result$y)),
    y     = rep(interp_result$y, each  = length(interp_result$x)),
    value = as.vector(interp_result$z))

  ggplot(df_grid, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    scale_fill_gradientn(
      colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
      limits = limits
    ) +
    coord_equal(expand = F, clip = "on") +
    ggtitle(title) +
    labs(x = "Easting", y = "Northing") +
    theme_paper(leg = "right")
}

plot_surface_mba <- function(mat, fixed_crd,
                             title        = NULL,
                             component    = 1,
                             grid_res     = 500,
                             limits       = NULL,
                             legend_label = NULL) {
  
  xyz <- data.frame(x = fixed_crd[, 1],
                    y = fixed_crd[, 2],
                    z = mat[, component])
  
  fit  <- MBA::mba.surf(xyz, no.X = grid_res, no.Y = grid_res, extend = TRUE)
  grid <- expand.grid(x = fit$xyz.est$x, y = fit$xyz.est$y)
  grid$z <- as.vector(fit$xyz.est$z)
  
  ggplot(grid, aes(x, y, fill = z)) +
    geom_raster(na.rm = TRUE) +
    scale_fill_gradientn(
      colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
      limits  = limits,
      oob     = scales::squish,
      guide   = guide_colorbar(
        title          = legend_label,
        title.position = "top",
        title.hjust    = 0.5
      )
    ) +
    coord_fixed(ratio = 1, clip = "on") + 
    ggtitle(title) +
    labs(x = "Easting", y = "Northing") +
    theme_paper(leg = "right") +
    theme(legend.title = element_text())
}
