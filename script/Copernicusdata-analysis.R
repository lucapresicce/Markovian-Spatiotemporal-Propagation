################################################################################
# COPERNICUS DATA ANALYSIS
################################################################################

# 0. packages -------------------------------------------------------------

library(spFFBS); library(spBPS); library(mniw)
library(fields);  library(dlm);  library(scoringRules)
library(nimble);  library(ggplot2); library(dplyr); library(tidyr)
library(gstat); library(spacetime)

has_inla <- requireNamespace("INLA", quietly = TRUE)
if (has_inla) library(INLA)


# 1. data -----------------------------------------------------------------

load("data/copernicus_data.RData")
load("data/copernicus_predictors.RData")

N_full <- dim(spatial_data)[1]
q      <- dim(spatial_data)[2]
t_full <- dim(spatial_data)[3]
p_raw  <- dim(spatial_predictors)[2]

tt <- 144
set.seed(42)
idx <- sample(1:N_full, 600, FALSE)

spatial_df   <- spatial_data[idx, , (t_full - tt + 1):t_full]
spatial_pred <- spatial_predictors[idx, , (t_full - tt + 1):t_full]
coords_all   <- coords[idx, ]

N <- dim(spatial_df)[1]
q <- dim(spatial_df)[2]
t <- dim(spatial_df)[3]

dates_vector <- seq(as.Date("1850-01-01"), by = "1 month", length.out = 1980)
dates_vector <- dates_vector[(1980 - tt + 1):1980]

design_matrix <- array(0, c(N, p_raw + 11, t))
for (i in 1:t) {
  P_t  <- scale(spatial_pred[, , i])
  L_t  <- matrix(0, N, 11)
  mese <- as.numeric(format(dates_vector[i], "%m"))
  if (mese != 12) L_t[, mese] <- 1
  design_matrix[, , i] <- cbind(P_t, L_t)
}
p <- dim(design_matrix)[2]

u <- 100; h <- 24; n <- N - u
set.seed(42)
train_idx   <- sample(1:N, n, FALSE)
test_idx    <- setdiff(1:N, train_idx)
tmax        <- t - h
forecast_id <- (tmax + 1):t

Ys   <- spatial_df[train_idx, , 1:tmax]
Xs   <- design_matrix[train_idx, , 1:tmax]
crds <- coords_all[train_idx, ]
Ds   <- spBPS:::arma_dist(crds)
D_mc <- as.matrix(dist(crds))

Yu   <- spatial_df[test_idx, , ]
Xu   <- design_matrix[test_idx, , ]
crdu <- coords_all[test_idx, ]

m0    <- matrix(0, n + p, q)
C0    <- rbind(cbind(diag(0.005, p), matrix(0, p, n)),
               cbind(matrix(0, n, p), exp(-Ds)))
prior <- list(m = m0, C = C0, nu = 3, Psi = diag(q))
L     <- 200

y_true_fcst <- spatial_df[train_idx, , forecast_id]   # n x q x h

cat(sprintf("n=%d | q=%d | tmax=%d | h=%d | u=%d | p=%d\n", n, q, tmax, h, u, p))


# 2. MNIW-DLM functions ---------------------------------------------------

lgamma_mq <- function(a, q) sum(lgamma(a - (seq_len(q) - 1) / 2))

mniw_kf <- function(tau, phi, t_end = tmax, store = FALSE) {
  a_v   <- 1 / tau - 1
  V     <- a_v * diag(n)                        # diagonal - Eq.(8) V_t(alpha)
  W_phi <- rbind(                               # state evolution noise
    cbind(diag(p),         matrix(0, p, n)),
    cbind(matrix(0, n, p), exp(-phi * D_mc))
  )
  C0_p <- rbind(
    cbind(diag(0.005, p),  matrix(0, p, n)),
    cbind(matrix(0, n, p), exp(-phi * D_mc))
  )
  m <- prior$m; C <- C0_p
  nu <- max(prior$nu, q + 2); Psi <- prior$Psi
  log_ml <- 0
  if (store) m_lst <- C_lst <- a_lst <- R_lst <- vector("list", t_end)
  for (tt_i in 1:t_end) {
    Pt  <- cbind(Xs[, , tt_i], diag(n)); Yt <- Ys[, , tt_i]
    a_t <- m
    R_t <- C + W_phi                            # prediction: add W at every step
    f_t <- Pt %*% a_t; Q_t <- Pt %*% R_t %*% t(Pt) + V
    Qc  <- tryCatch(chol(Q_t + diag(n) * 1e-10), error = function(e) diag(n))
    Qi  <- chol2inv(Qc); e_t <- Yt - f_t
    Psi_n <- Psi + t(e_t) %*% Qi %*% e_t; nu_n <- nu + n
    lQ  <- 2 * sum(log(diag(Qc)))
    lP0 <- as.numeric(determinant(Psi,   log = TRUE)$modulus)
    lPn <- as.numeric(determinant(Psi_n, log = TRUE)$modulus)
    log_ml <- log_ml + (-0.5 * n * q * log(pi) +
                          lgamma_mq(nu_n / 2, q) - lgamma_mq(nu / 2, q) +
                          0.5 * nu * lP0 - 0.5 * nu_n * lPn - 0.5 * q * lQ)
    Kt <- R_t %*% t(Pt) %*% Qi
    m  <- a_t + Kt %*% e_t
    C  <- (diag(n + p) - Kt %*% Pt) %*% R_t; C <- (C + t(C)) / 2
    nu <- nu_n; Psi <- Psi_n
    if (store) {
      a_lst[[tt_i]] <- a_t; R_lst[[tt_i]] <- R_t
      m_lst[[tt_i]] <- m;   C_lst[[tt_i]] <- C
    }
  }
  out <- list(m = m, C = C, nu = nu, Psi = Psi, log_ml = log_ml)
  if (store) out <- c(out, list(m_lst = m_lst, C_lst = C_lst,
                                a_lst = a_lst, R_lst = R_lst))
  out
}

mniw_kf_step <- function(m_p, C_p, nu_p, Psi_p, tt_i, tau, phi) {
  a_v   <- 1 / tau - 1
  V     <- a_v * diag(n)
  W_phi <- rbind(
    cbind(diag(p),         matrix(0, p, n)),
    cbind(matrix(0, n, p), exp(-phi * D_mc))
  )
  Pt  <- cbind(Xs[, , tt_i], diag(n)); Yt <- Ys[, , tt_i]
  a_t <- m_p
  R_t <- C_p + W_phi
  f_t <- Pt %*% a_t; Q_t <- Pt %*% R_t %*% t(Pt) + V
  tryCatch({
    Qc    <- chol(Q_t + diag(n) * 1e-10); Qi <- chol2inv(Qc); e_t <- Yt - f_t
    Psi_n <- Psi_p + t(e_t) %*% Qi %*% e_t; nu_n <- nu_p + n
    lQ    <- 2 * sum(log(diag(Qc)))
    lP0   <- as.numeric(determinant(Psi_p,  log = TRUE)$modulus)
    lPn   <- as.numeric(determinant(Psi_n,  log = TRUE)$modulus)
    lp    <- -0.5 * n * q * log(pi) +
      lgamma_mq(nu_n / 2, q) - lgamma_mq(nu_p / 2, q) +
      0.5 * nu_p * lP0 - 0.5 * nu_n * lPn - 0.5 * q * lQ
    Kt  <- R_t %*% t(Pt) %*% Qi
    m_n <- a_t + Kt %*% e_t
    C_n <- (diag(n + p) - Kt %*% Pt) %*% R_t; C_n <- (C_n + t(C_n)) / 2
    list(m = m_n, C = C_n, nu = nu_n, Psi = Psi_n, lp = lp, ok = TRUE)
  }, error = function(e)
    list(m = m_p, C = C_p, nu = nu_p, Psi = Psi_p, lp = -Inf, ok = FALSE))
}

mniw_forecast_samp <- function(kf, tau, phi, L_ = L) {
  arr   <- array(0, c(n, q, L_, h))
  a_v   <- 1 / tau - 1
  W_phi <- rbind(
    cbind(diag(p),         matrix(0, p, n)),
    cbind(matrix(0, n, p), exp(-phi * D_mc))
  )
  for (l in 1:L_) {
    Sig_l <- mniw::riwish(1, Psi = kf$Psi, nu = kf$nu)
    LcS   <- t(chol(Sig_l + diag(q) * 1e-8))
    for (hh in 1:h) {
      R_h   <- kf$C + hh * W_phi
      LcR_h <- tryCatch(t(chol(R_h + diag(n + p) * 1e-8)),
                        error = function(e) diag(n + p))
      Theta_h <- kf$m + LcR_h %*% matrix(rnorm((n + p) * q), n + p, q) %*% LcS
      Pt_h    <- cbind(design_matrix[train_idx, , tmax + hh], diag(n))
      arr[, , l, hh] <- Pt_h %*% Theta_h +
        sqrt(a_v) * matrix(rnorm(n * q), n, q) %*% LcS
    }
  }
  arr
}

mniw_krig_samp <- function(kf, tau, phi, t_pred = tmax, L_ = L) {
  Df  <- as.matrix(dist(rbind(crds, crdu)))
  Kf  <- exp(-phi * Df)
  K11 <- Kf[1:n, 1:n];                    diag(K11) <- diag(K11) + 1e-6
  K12 <- Kf[1:n, (n + 1):(n + u)]
  K22 <- Kf[(n + 1):(n + u), (n + 1):(n + u)]; diag(K22) <- diag(K22) + 1e-6
  Ak  <- t(K12) %*% chol2inv(chol(K11))
  Ck  <- K22 - Ak %*% K12
  LcK <- tryCatch(t(chol(Ck + diag(u) * 1e-8)), error = function(e) diag(u))
  W_phi <- rbind(
    cbind(diag(p),         matrix(0, p, n)),
    cbind(matrix(0, n, p), exp(-phi * D_mc))
  )
  steps_ahead <- max(t_pred - tmax, 0L)
  C_pred <- kf$C + steps_ahead * W_phi
  LcW    <- tryCatch(t(chol(C_pred + diag(n + p) * 1e-8)),
                     error = function(e) diag(n + p))
  out <- array(0, c(u, q, L_))
  for (l in 1:L_) {
    Sig_l   <- mniw::riwish(1, Psi = kf$Psi, nu = kf$nu)
    LcS     <- t(chol(Sig_l + diag(q) * 1e-8))
    W_T     <- kf$m + LcW %*% matrix(rnorm((n + p) * q), n + p, q) %*% LcS
    beta_T  <- W_T[1:p, , drop = FALSE]
    omega_T <- W_T[(p + 1):(p + n), , drop = FALSE]
    omega_u <- Ak %*% omega_T
    X_te    <- design_matrix[test_idx, , t_pred]
    out[, , l] <- X_te %*% beta_T + omega_u +
      LcK %*% matrix(rnorm(u * q), u, q) %*% LcS
  }
  out
}


# 3. evaluation setup -----------------------------------------------------

# In-sample: tmax; out-of-sample: h = 1, 6, 12, 24

sp_is  <- tmax
sp_os  <- c(1, 6, 12, 24)
T_sp   <- 1 + length(sp_os)

sp_labels <- c(paste0("IS_t", sp_is), paste0("OS_h", sp_os))

y_true_sp_all <- array(0, c(u, q, T_sp))
y_true_sp_all[, , 1] <- spatial_df[test_idx, , tmax]
for (k in seq_along(sp_os)) {
  y_true_sp_all[, , 1 + k] <- spatial_df[test_idx, , tmax + sp_os[k]]
  }

sp_all_fixed <- function(tau, phi, kf_tm, L_ = L) {
  arr <- array(0, c(u, q, L_, T_sp))
  arr[, , , 1] <- mniw_krig_samp(kf_tm, tau, phi, tmax,            L_)
  for (k in seq_along(sp_os))
    arr[, , , 1 + k] <- mniw_krig_samp(kf_tm, tau, phi, tmax + sp_os[k], L_)
  arr
}


# 4. metric functions -----------------------------------------------------

rmspe_fn <- function(yt, yp) {sqrt(mean((as.numeric(yt) - as.numeric(yp))^2))}

cov95_fn <- function(yt, yl, yu) {
  mean(as.numeric(yt) >= as.numeric(yl) & as.numeric(yt) <= as.numeric(yu))
  }

is_fn <- function(yt, yl, yu, a = 0.05) {
  yt <- as.numeric(yt); yl <- as.numeric(yl); yu <- as.numeric(yu)
  mean(yu - yl + (2 / a) * pmax(yl - yt, 0) + (2 / a) * pmax(yt - yu, 0))
}

es_fn <- function(yt, samp) {
  d1 <- dim(yt)[1]; d3 <- dim(yt)[3]
  vals <- numeric(d1 * d3); k <- 0
  for (t_ in 1:d3) for (s in 1:d1) {
    k <- k + 1
    vals[k] <- scoringRules::es_sample(y = yt[s, , t_], dat = samp[s, , , t_])
  }
  mean(vals)
}

compute_metrics <- function(yt, samp) {
  ym <- apply(samp, c(1, 2, 4), mean)
  yl <- apply(samp, c(1, 2, 4), quantile, 0.025)
  yu <- apply(samp, c(1, 2, 4), quantile, 0.975)
  list(rmspe  = rmspe_fn(yt, ym),
       cov95  = cov95_fn(yt, yl, yu),
       is     = is_fn(yt, yl, yu),
       es     = es_fn(yt, samp))
}

metrics_by_h <- function(tag, yt, samp) {
  h_ <- dim(yt)[3]
  do.call(rbind, lapply(1:h_, function(hh) {
    yh <- yt[, , hh]; sh <- samp[, , , hh]
    ym <- apply(sh, c(1, 2), mean)
    yl <- apply(sh, c(1, 2), quantile, 0.025)
    yu <- apply(sh, c(1, 2), quantile, 0.975)
    y3 <- array(yh, c(dim(yh), 1)); s4 <- array(sh, c(dim(sh), 1))
    data.frame(method = tag, h = hh,
               rmspe  = rmspe_fn(yh, ym), cov95 = cov95_fn(yh, yl, yu),
               is     = is_fn(yh, yl, yu),
               es     = es_fn(y3, s4))
  }))
}

metrics_by_sp <- function(tag, yt, samp) {
  do.call(rbind, lapply(1:T_sp, function(k) {
    yk <- yt[, , k]; sk <- samp[, , , k]
    ym <- apply(sk, c(1, 2), mean)
    yl <- apply(sk, c(1, 2), quantile, 0.025)
    yu <- apply(sk, c(1, 2), quantile, 0.975)
    y3 <- array(yk, c(dim(yk), 1)); s4 <- array(sk, c(dim(sk), 1))
    data.frame(method = tag, time = sp_labels[k],
               rmspe  = rmspe_fn(yk, ym), cov95 = cov95_fn(yk, yl, yu),
               is     = is_fn(yk, yl, yu),
               es     = es_fn(y3, s4))
  }))
}

print_met <- function(tag, mT, mS, secs) {
  cat(sprintf("  %s [T] RMSPE=%.4f COV95=%.4f IS=%.4f ES=%.4f\n",
              tag, mT$rmspe, mT$cov95, mT$is, mT$es))
  cat(sprintf("  %s [S] RMSPE=%.4f COV95=%.4f IS=%.4f ES=%.4f | %.0fs\n",
              tag, mS$rmspe, mS$cov95, mS$is, mS$es, secs))
}


# 5. DYNBPS ---------------------------------------------------------------

Rcpp::sourceCpp("script/FFBS-DYNBPS-struct-v2.cpp")
source("script/FFBS-DYNBPS-v2.R")

# load variogram-informed grid
load("data/variogram_informed_grid.Rdata")
par_grid <- spBPS:::expand_grid_cpp(rev(tau_seq), rev(phi_seq))
J        <- nrow(par_grid)

sp_pred_idx <- c(sp_is, sp_is+sp_os)

cat(">>> [1/5] DYNBPS...\n"); t0 <- proc.time(); set.seed(42)
out_DYN <- spFFBS_v2(
  Y = Ys, X = Xs, X_all = design_matrix[train_idx, , ], D = Ds,
  grid    = list(tau = tau_seq, phi = phi_seq),
  prior   = prior, L = L,
  do_BS   = TRUE, do_forecast = TRUE, do_spatial = TRUE, tnew = h,
  spatial = list(crd = crds, crdtilde = crdu, Xtilde = Xu[, , sp_pred_idx], t = sp_pred_idx),
  n_threads = 5L)
time_DYN <- (proc.time() - t0)["elapsed"]

Y_fcst_DYN   <- out_DYN$forecast$Y_pred[, , , -(1:tmax), drop = FALSE]
W_global_DYN <- out_DYN$Wglobal
Ysp_DYN_all <- simplify2array(lapply(out_DYN$spatial, `[[`, "Y"))

met_DYN_T <- compute_metrics(y_true_fcst, Y_fcst_DYN)
met_DYN_S <- compute_metrics(y_true_sp_all, Ysp_DYN_all)
ph_DYN    <- metrics_by_h( "DYNBPS", y_true_fcst, Y_fcst_DYN)
ps_DYN    <- metrics_by_sp("DYNBPS", y_true_sp_all, Ysp_DYN_all)
print_met("DYNBPS", met_DYN_T, met_DYN_S, time_DYN)


# 6. EB -------------------------------------------------------------------

cat("\n>>> [2/5] EB-static...\n"); t0 <- proc.time()

neg_lml_EB <- function(par) {
  tau <- plogis(par[1]) * 0.9999 + 1e-5; phi <- exp(par[2])
  -mniw_kf(tau, phi, tmax)$log_ml
}
opt_EB <- optim(c(0, log(1)), neg_lml_EB,
                method = "L-BFGS-B",
                lower  = c(qlogis(0.001), log(0.001)),
                upper  = c(qlogis(0.9999), log(100.0)),
                control = list(maxit = 500))
tau_EB <- plogis(opt_EB$par[1]) * 0.9999 + 1e-5
phi_EB <- exp(opt_EB$par[2])
cat(sprintf("  EB: tau=%.4f  phi=%.4f\n", tau_EB, phi_EB))

set.seed(42)
kf_EB      <- mniw_kf(tau_EB, phi_EB, tmax)
Y_fcst_EB  <- mniw_forecast_samp(kf_EB, tau_EB, phi_EB, L)
Ysp_EB_all <- sp_all_fixed(tau_EB, phi_EB, kf_EB, L)
time_EB    <- (proc.time() - t0)["elapsed"]

met_EB_T <- compute_metrics(y_true_fcst, Y_fcst_EB)
met_EB_S <- compute_metrics(y_true_sp_all, Ysp_EB_all)
ph_EB    <- metrics_by_h( "EB-static", y_true_fcst, Y_fcst_EB)
ps_EB    <- metrics_by_sp("EB-static", y_true_sp_all, Ysp_EB_all)
print_met("EB-st", met_EB_T, met_EB_S, time_EB)


# 7. MCMC -----------------------------------------------------------------

cat("\n>>> [3/5] MCMC-nimble...\n")

mniw_ml_R_wrap <- function(tau, phi)
  tryCatch(mniw_kf(tau, phi, tmax)$log_ml, error = function(e) -1e6)

mniw_ml_nf <- nimbleRcall(
  prototype  = function(tau = double(0), phi = double(0)) {},
  returnType = double(0), Rfun = "mniw_ml_R_wrap")

code_nm <- nimbleCode({
  tau        ~ dunif(0.001, 0.9999)
  phi        ~ T(dgamma(2, 5), 0.01, 10)
  log_ml_val <- mniw_ml_nf(tau, phi)
  zeros[1]   ~ dpois(-log_ml_val)
})
model_nm  <- nimbleModel(code_nm, data = list(zeros = 0L),
                         inits = list(tau = tau_EB, phi = max(phi_EB, 0.05)),
                         calculate = FALSE)
cmodel_nm <- compileNimble(model_nm)
conf_nm   <- configureMCMC(cmodel_nm, nodes = character(0))
conf_nm$addSampler("tau", type = "RW",
                   control = list(adaptive = TRUE, adaptInterval = 100, scale = 0.15))
conf_nm$addSampler("phi", type = "RW",
                   control = list(adaptive = TRUE, adaptInterval = 100, scale = 0.20, log = TRUE))
mcmc_nm  <- buildMCMC(conf_nm)
cmcmc_nm <- compileNimble(mcmc_nm, project = cmodel_nm, resetFunctions = TRUE)

t0_nm    <- proc.time()
samps_nm <- runMCMC(cmcmc_nm, niter = 2000, nburnin = 1000,
                    progressBar = TRUE, summary = FALSE)
cat(sprintf("\n  MCMC: tau=%.4f  phi=%.4f\n",
            mean(samps_nm[, "tau"]), mean(samps_nm[, "phi"])))

keep_nm    <- sample(nrow(samps_nm), L)
set.seed(42)
Y_fcst_NM  <- array(0, c(n, q, L, h))
Ysp_NM_all <- array(0, c(u, q, L, T_sp))
for (l in seq_along(keep_nm)) {
  tau_l <- samps_nm[keep_nm[l], "tau"]; phi_l <- samps_nm[keep_nm[l], "phi"]
  kf_l  <- mniw_kf(tau_l, phi_l, tmax)
  Y_fcst_NM[, , l, ]    <- mniw_forecast_samp(kf_l, tau_l, phi_l, 1)[, , 1, ]
  Ysp_NM_all[, , l, 1]  <- mniw_krig_samp(kf_l, tau_l, phi_l, tmax, 1)[, , 1]
  for (k in seq_along(sp_os))
    Ysp_NM_all[, , l, 1 + k] <- mniw_krig_samp(kf_l, tau_l, phi_l, tmax + sp_os[k], 1)[, , 1]
}
time_NM <- (proc.time() - t0_nm)["elapsed"]

met_NM_T <- compute_metrics(y_true_fcst, Y_fcst_NM)
met_NM_S <- compute_metrics(y_true_sp_all, Ysp_NM_all)
ph_NM    <- metrics_by_h( "MCMC", y_true_fcst, Y_fcst_NM)
ps_NM    <- metrics_by_sp("MCMC", y_true_sp_all, Ysp_NM_all)
print_met("MCMC", met_NM_T, met_NM_S, time_NM); gc()


# 8. EOF ------------------------------------------------------------------

cat("\n>>> [5a] EOF+DLM...\n"); t0_eof <- proc.time()
K_eof <- 10

Ymat_full  <- matrix(0, N * q, tmax)
for (i in 1:tmax) Ymat_full[, i] <- as.numeric(spatial_df[, , i])
Ymean_full <- rowMeans(Ymat_full)
sv         <- svd(Ymat_full - Ymean_full, nu = K_eof, nv = K_eof)
U_eof_full <- sv$u                              # (N*q) x K
sc_eof     <- diag(sv$d[1:K_eof]) %*% t(sv$v)  # K x tmax

cat(sprintf("  EOF variance explained: %.1f%%\n",
            100 * sum(sv$d[1:K_eof]^2) / sum(sv$d^2)))

train_rows_eof <- as.integer(outer(train_idx, (0:(q-1)) * N, "+"))
test_rows_eof  <- as.integer(outer(test_idx,  (0:(q-1)) * N, "+"))

U_tr     <- U_eof_full[train_rows_eof, ]   # (n*q) x K
U_te     <- U_eof_full[test_rows_eof,  ]   # (u*q) x K
Ymean_tr <- Ymean_full[train_rows_eof]     # n*q
Ymean_te <- Ymean_full[test_rows_eof]      # u*q

U_krow <- matrix(0, K_eof, n)
for (k in 1:K_eof) U_krow[k, ] <- rowSums(matrix(U_tr[, k], n, q))

U_krow_test <- matrix(0, K_eof, u)
for (k in 1:K_eof) U_krow_test[k, ] <- rowSums(matrix(U_te[, k], u, q))

Xproj_tr <- array(0, c(K_eof, p, tmax))
for (tt_i in 1:tmax) Xproj_tr[, , tt_i] <- U_krow %*% Xs[, , tt_i]

Xproj_fc <- array(0, c(K_eof, p, h))
for (hh in 1:h) Xproj_fc[, , hh] <- U_krow %*% design_matrix[train_idx, , tmax + hh]

mniw_kf_eof <- function(tau, t_end = tmax) {
  a_v <- 1 / tau - 1
  V_K <- diag(K_eof) * a_v + diag(K_eof) * 1e-6
  C0K <- rbind(cbind(diag(0.005, p),    matrix(0, p, K_eof)),
               cbind(matrix(0, K_eof, p), diag(K_eof)))
  m <- matrix(0, K_eof + p, 1); C <- C0K
  nu <- 3L; Psi <- matrix(1); log_ml <- 0
  for (tt_i in 1:t_end) {
    Pt  <- cbind(Xproj_tr[, , tt_i], diag(K_eof))
    Yt  <- matrix(sc_eof[, tt_i], K_eof, 1)
    a_t <- m; R_t <- C
    f_t <- Pt %*% a_t; Q_t <- Pt %*% R_t %*% t(Pt) + V_K
    Qc  <- tryCatch(chol(Q_t + diag(K_eof) * 1e-10), error = function(e) diag(K_eof))
    Qi  <- chol2inv(Qc); e_t <- Yt - f_t
    Psi_n <- as.matrix(Psi + t(e_t) %*% Qi %*% e_t); nu_n <- nu + K_eof
    lQ    <- 2 * sum(log(diag(Qc)))
    log_ml <- log_ml + (-0.5 * K_eof * log(pi) +
                          lgamma(nu_n / 2) - lgamma(nu / 2) +
                          0.5 * nu   * log(max(as.numeric(Psi),   1e-300)) -
                          0.5 * nu_n * log(max(as.numeric(Psi_n), 1e-300)) -
                          0.5 * lQ)
    Kt <- R_t %*% t(Pt) %*% Qi
    m  <- a_t + Kt %*% e_t
    C  <- (diag(K_eof + p) - Kt %*% Pt) %*% R_t; C <- (C + t(C)) / 2
    nu <- nu_n; Psi <- Psi_n
  }
  list(m = m, C = C, nu = nu, Psi = Psi, log_ml = log_ml)
}

opt_eof <- optim(0, function(par) { tau <- plogis(par) * 0.9999 + 1e-5
-mniw_kf_eof(tau)$log_ml }, method = "Brent", lower = -6, upper = 6)
tau_eof <- plogis(opt_eof$par) * 0.9999 + 1e-5
kf_eof  <- mniw_kf_eof(tau_eof, tmax)
cat(sprintf("  EOF: tau=%.4f\n", tau_eof))

LcC_eof <- tryCatch(t(chol(kf_eof$C + diag(K_eof + p) * 1e-8)),
                    error = function(e) diag(K_eof + p))

set.seed(42)
Y_fcst_EOF  <- array(0, c(n, q, L, h))
Ysp_EOF_all <- array(0, c(u, q, L, T_sp))

for (l in 1:L) {
  sig2_l <- 1 / rgamma(1, shape = kf_eof$nu / 2,
                       rate  = as.numeric(kf_eof$Psi) / 2)
  sd_l   <- sqrt(sig2_l)
  W_T_l  <- kf_eof$m + sd_l * (LcC_eof %*% rnorm(K_eof + p))
  B_K_l  <- W_T_l[1:p]
  xi_T_l <- W_T_l[(p + 1):(p + K_eof)]
  
  # Temporal forecast: reconstruct via training rows of U
  for (hh in 1:h) {
    s_hat <- as.numeric(Xproj_fc[, , hh] %*% B_K_l) + xi_T_l
    Y_fcst_EOF[, , l, hh] <- matrix(as.numeric(U_tr %*% s_hat) + Ymean_tr, n, q)
  }
  
  # Spatial: reconstruct via test rows of U - no external interpolation
  for (k_ in 1:T_sp) {
    t_k <- if (k_ == 1) tmax else tmax + sp_os[k_ - 1]
    Xp  <- U_krow_test %*% design_matrix[test_idx, , t_k]   # K x p
    s_k <- as.numeric(Xp %*% B_K_l) + xi_T_l
    Ysp_EOF_all[, , l, k_] <- matrix(as.numeric(U_te %*% s_k) + Ymean_te, u, q) +
      matrix(rnorm(u * q, 0, sd_l), u, q)
  }
}
time_EOF <- (proc.time() - t0_eof)["elapsed"]

met_EOF_T <- compute_metrics(y_true_fcst, Y_fcst_EOF)
met_EOF_S <- compute_metrics(y_true_sp_all, Ysp_EOF_all)
ph_EOF    <- metrics_by_h( "EOF+DLM", y_true_fcst, Y_fcst_EOF)
ps_EOF    <- metrics_by_sp("EOF+DLM", y_true_sp_all, Ysp_EOF_all)
print_met("EOF", met_EOF_T, met_EOF_S, time_EOF)


# 9. INLA -----------------------------------------------------------------

cat("\n>>> [5b] INLA-multiv...\n")

if (!has_inla) {
  met_IM_T <- met_IM_S <- ph_IM <- ps_IM <- NULL; time_IM <- NA
  cat("  INLA not available.\n")
} else {
  t0_im  <- proc.time(); L_inla <- 50
  mesh_u <- INLA::inla.mesh.2d(loc = crds, max.edge = c(20, 40), cutoff = 5)
  spde_u <- INLA::inla.spde2.pcmatern(mesh_u,
                                      prior.range = c(0.1, 0.5),
                                      prior.sigma = c(1.0, 0.5))
  n_spde <- spde_u$n.spde
  cat(sprintf("  Mesh nodes: %d\n", n_spde))
  
  n_obs_m <- n * q * tmax
  y_m  <- numeric(n_obs_m); sp_m <- gr_m <- rp_m <- integer(n_obs_m)
  X_full <- matrix(0, n_obs_m, p * q)
  row_idx <- 0
  for (v in 1:q) for (tt_i in 1:tmax) {
    rows <- row_idx + seq_len(n)
    y_m[rows]  <- Ys[, v, tt_i]; sp_m[rows] <- seq_len(n)
    gr_m[rows] <- tt_i;          rp_m[rows] <- v
    X_full[rows, ((v - 1) * p + 1):(v * p)] <- Xs[, , tt_i]
    row_idx <- row_idx + n
  }
  A_m    <- INLA::inla.spde.make.A(mesh_u, loc = crds[sp_m, , drop = FALSE],
                                   group = gr_m, n.group = tmax,
                                   repl = rp_m, n.repl = q)
  iset_m <- INLA::inla.spde.make.index("space", n.spde = n_spde,
                                       n.group = tmax, n.repl = q)
  stk_m  <- INLA::inla.stack(data    = list(y = y_m),
                             A       = list(A_m, 1),
                             effects = list(iset_m, list(X_full = X_full)),
                             tag = "obs")
  
  res_m <- tryCatch(
    INLA::inla(
      y ~ -1 + X_full + f(space, model = spde_u, group = space.group,
                          control.group = list(model = "ar1"),
                          replicate = space.repl),
      data = INLA::inla.stack.data(stk_m), family = "gaussian",
      control.predictor = list(A = INLA::inla.stack.A(stk_m), compute = TRUE),
      control.compute  = list(config = TRUE), verbose = FALSE),
    error = function(e) { cat("  INLA failed:", e$message, "\n"); NULL })
  
  if (is.null(res_m)) {
    met_IM_T <- met_IM_S <- ph_IM <- ps_IM <- NULL
    time_IM  <- (proc.time() - t0_im)["elapsed"]
  } else {
    A_tr    <- INLA::inla.spde.make.A(mesh_u, loc = crds)
    A_te    <- INLA::inla.spde.make.A(mesh_u, loc = crdu)
    samps_m <- INLA::inla.posterior.sample(L_inla, res_m)
    mvar_m  <- tryCatch(
      INLA::inla.spde.result(res_m, "space", spde_u)$summary.log.variance.nominal["mean"],
      error = function(e) 0)
    sf_m <- exp(as.numeric(mvar_m) / 2)
    
    Y_fcst_IM  <- array(0, c(n, q, L_inla, h))
    Ysp_IM_all <- array(0, c(u, q, L_inla, T_sp))
    
    for (s in 1:L_inla) {
      lat      <- samps_m[[s]]$latent
      rho      <- as.numeric(samps_m[[s]]$hyperpar["GroupRho for space"])
      sig2     <- 1 / as.numeric(samps_m[[s]]$hyperpar[
        "Precision for the Gaussian observations"])
      inn_sd   <- sf_m * sqrt(max(1 - rho^2, 0))
      beta_all <- lat[grep("^X_full", rownames(lat)), 1]
      sp_idx   <- which(grepl("^space:", rownames(lat)))
      
      for (v in 1:q) {
        beta_v <- beta_all[((v - 1) * p + 1):(v * p)]
        off_v  <- (v - 1) * n_spde * tmax
        sp_Tv  <- lat[sp_idx[(off_v + (tmax - 1) * n_spde + 1):(off_v + tmax * n_spde)], 1]
        
        sp_cur <- sp_Tv
        for (hh in 1:h) {
          sp_cur <- rho * sp_cur + rnorm(n_spde, 0, inn_sd)
          mu_h   <- as.numeric(A_tr %*% sp_cur) +
            as.numeric(design_matrix[train_idx, , tmax + hh] %*% beta_v)
          Y_fcst_IM[, v, s, hh] <- rnorm(n, mu_h, sqrt(sig2))
        }
        mu_sp <- as.numeric(A_te %*% sp_Tv) +
          as.numeric(design_matrix[test_idx, , tmax] %*% beta_v)
        Ysp_IM_all[, v, s, 1] <- rnorm(u, mu_sp, sqrt(sig2))
        for (k in seq_along(sp_os)) {
          sp_hh <- sp_Tv
          for (step in 1:sp_os[k]) sp_hh <- rho * sp_hh + rnorm(n_spde, 0, inn_sd)
          mu_k  <- as.numeric(A_te %*% sp_hh) +
            as.numeric(design_matrix[test_idx, , tmax + sp_os[k]] %*% beta_v)
          Ysp_IM_all[, v, s, 1 + k] <- rnorm(u, mu_k, sqrt(sig2))
        }
      }
    }
    time_IM  <- (proc.time() - t0_im)["elapsed"]
    met_IM_T <- compute_metrics(y_true_fcst, Y_fcst_IM)
    met_IM_S <- compute_metrics(y_true_sp_all, Ysp_IM_all)
    ph_IM    <- metrics_by_h( "INLA-multiv", y_true_fcst, Y_fcst_IM)
    ps_IM    <- metrics_by_sp("INLA-multiv", y_true_sp_all, Ysp_IM_all)
    print_met("INLA-m", met_IM_T, met_IM_S, time_IM)
  }; gc()
}


# 10. tables --------------------------------------------------------------

method_names <- c("DYNBPS", "EB-static", "MCMC", "EOF+DLM", "INLA-multiv")
mets_T  <- list(met_DYN_T, met_EB_T, met_NM_T, met_EOF_T, met_IM_T)
mets_S  <- list(met_DYN_S, met_EB_S, met_NM_S, met_EOF_S, met_IM_S)
times_s <- c(time_DYN, time_EB, time_NM, time_EOF, time_IM)

gf <- function(m, f) if (!is.null(m)) round(m[[f]], 5) else NA

make_df <- function(mets)
  data.frame(Method   = method_names,
             RMSPE    = sapply(mets, gf, "rmspe"),
             Coverage = sapply(mets, gf, "cov95"),
             IS       = sapply(mets, gf, "is"),
             ES       = sapply(mets, gf, "es"),
             Time_s   = round(times_s, 1))

tab_T <- make_df(mets_T); tab_S <- make_df(mets_S)
cat("\n===== Table A: Temporal =====\n"); print(tab_T, row.names = FALSE)
cat("\n===== Table B: Spatial =====\n");  print(tab_S, row.names = FALSE)


# 11. plots ---------------------------------------------------------------

library(MBA); library(sf); library(rnaturalearth)
library(patchwork); library(ggh4x); library(purrr)
library(RColorBrewer); library(scales); library(tibble)


# 11.0.  Shared helpers ---------------------------------------------------

theme_paper <- function(base = 14, leg = "none") {
  theme_bw(base_size = base) %+replace% theme(
    title            = element_text(size = base + 2, face = "bold"),
    strip.text       = element_text(size = base,     face = "bold"),
    strip.background = element_rect(fill = "#eceff1", colour = "grey70"),
    axis.text.x      = element_text(size = base - 2, hjust = 0.5),
    axis.text.y      = element_text(size = base - 2, hjust = 0.5),
    axis.title       = element_text(size = base,     face = "bold"),
    axis.title.x     = element_text(size = base,     face = "bold"),
    axis.title.y     = element_text(size = base,     face = "bold", angle = 90),
    legend.text      = element_text(size = base + 2, face = "bold"),
    legend.title     = element_blank(),
    legend.key.size  = unit(1.3, "lines"),
    legend.position  = leg,
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )
}

europe_sf    <- ne_countries(continent = "Europe", scale = "medium", returnclass = "sf")
var_labels   <- c("Temp", "Rain", "Wind", "Evps")
var_palettes <- c("RdBu", "C", "E", "D")   # RdBu for Temp; viridis C/E/D for rest

geo_surface <- function(lon, lat, z,
                        title = NULL, limits = NULL,
                        pal = "RdBu", legend_label = NULL,
                        grid_res = 300) {
  ok  <- complete.cases(data.frame(lon, lat, z))
  fit <- MBA::mba.surf(data.frame(x = lon[ok], y = lat[ok], z = z[ok]),
                       no.X = grid_res, no.Y = grid_res, extend = TRUE)
  grd <- expand.grid(Longitude = fit$xyz.est$x, Latitude = fit$xyz.est$y)
  grd$z <- as.vector(fit$xyz.est$z)
  
  fscale <- if (pal == "RdBu") {
    scale_fill_gradientn(
      colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
      limits = limits, oob = scales::squish, na.value = "transparent",
      guide  = guide_colorbar(title = legend_label,
                              title.position = "top", title.hjust = 0.5))
  } else {
    scale_fill_viridis_c(
      option = pal, direction = -1, limits = limits,
      oob = scales::squish, na.value = "transparent",
      guide = guide_colorbar(title = legend_label,
                             title.position = "top", title.hjust = 0.5))
  }
  
  ggplot() +
    geom_raster(data = grd,
                aes(x = Longitude, y = Latitude, fill = z), interpolate = TRUE) +
    geom_sf(data = europe_sf, fill = NA, colour = "black", linewidth = 0.2) +
    fscale +
    coord_sf(xlim = range(lon[ok]), ylim = range(lat[ok]), expand = FALSE) +
    ggtitle(title) +
    theme_paper(leg = "right") +
    theme(legend.title = element_text(face = "bold", size = 11),
          axis.title   = element_blank(),
          axis.text    = element_text(size = 9))
}


# 11.1.  Posterior summaries DYNBPS ---------------------------------------

Y_postmean <- apply(out_DYN$forecast$Y_pred, c(1, 2, 4), mean)
Y_upp      <- apply(out_DYN$forecast$Y_pred, c(1, 2, 4), quantile, 0.975)
Y_low      <- apply(out_DYN$forecast$Y_pred, c(1, 2, 4), quantile, 0.025)


# 11.2.  Figure 6a - Europe map with 4 forecast sites ---------------------

smp_loc <- c(419, 222, 67, 97)

lon_tr  <- crds[, 1]; lat_tr <- crds[, 2]
slabs   <- paste0("S", seq_along(smp_loc))

gg_points <- ggplot() +
  geom_sf(data   = europe_sf,
          fill   = rgb(190/255, 220/255, 180/255, 0.45),
          colour = "black", linewidth = 0.12) +
  geom_point(aes(x = lon_tr[smp_loc], y = lat_tr[smp_loc]),
             colour = "red", size = 16, pch = "*") +
  geom_text(aes(x     = lon_tr[smp_loc],
                y     = lat_tr[smp_loc],
                label = slabs),
            nudge_x = -2.2, nudge_y = 1.2,
            colour = "red", size = 5.5, fontface = "bold") +
  coord_sf(xlim = range(lon_tr) * c(0.97, 1.05),
           ylim = range(lat_tr) * c(0.97, 1.03),
           expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_paper(base = 18) +
  theme(panel.background = element_rect(fill = "#d6eaf8"),
        panel.border     = element_rect(colour = "grey40", linewidth = 0.4))


# 11.3.  Figure 6b - time-series forecast with 95% credible interval ------

truth_full <- spatial_df[train_idx, , ]   # [n, q, tmax+h]

make_tsdf <- function(i, j) {
  purrr::map_dfr(seq_len(q), function(v) {
    tibble(
      time     = seq_len(tmax + h),
      smp      = factor(slabs[i], levels = slabs),
      var      = v,
      truth    = truth_full[j, v, ],
      low      = Y_low[j, v, ],
      upp      = Y_upp[j, v, ],
      postmean = Y_postmean[j, v, ]
    )
  })
}
ts_df <- purrr::map2_dfr(seq_along(smp_loc), smp_loc, make_tsdf)

gg_lines <- ggplot(ts_df, aes(x = time)) +
  geom_ribbon(aes(ymin = low, ymax = upp, fill = factor(var)),
              alpha = 0.35, colour = "red", linetype = "dashed",
              linewidth = 0.12, show.legend = TRUE) +
  geom_line(aes(y = truth),    colour = "black", linewidth = 0.18) +
  geom_line(aes(y = postmean), colour = "red",   linewidth = 0.18) +
  geom_vline(xintercept = tmax + 0.5,
             linetype = "dotted", colour = "grey40", linewidth = 0.45) +
  facet_grid2(rows     = vars(smp), cols = vars(var),
              scales   = "free",  independent = "all",
              space    = "free",  axes = "all",
              labeller = labeller(var = setNames(var_labels, seq_len(q)))) +
  scale_fill_brewer(palette = "Accent", labels = var_labels) +
  scale_x_continuous(breaks = c(1, 30, 60, 90, 120, 144)) +
  labs(x = "Time (months)", y = NULL, fill = "variables") +
  guides(fill = guide_legend(override.aes = list(colour = NA, linetype = 0))) +
  theme_paper(leg = "bottom", base = 18) +
  theme(strip.text.y = element_text(angle = 0, size = 16),
        axis.title.y = element_blank(),
        legend.text  = element_text(size = 18, face = "bold"))

ggsave("plots/copernicus_temporal_forecast_points.png", gg_points, dpi = 300)
ggsave("plots/copernicus_temporal_forecast_lines.png", gg_lines, width = 14, height = 7, dpi = 300)
cat("Fig 6 saved.\n")


# 11.4.  Figure 7 ---------------------------------------------------------

geo_surface <- function(lon, lat, z,
                        title = NULL, limits = NULL,
                        pal = "RdBu", legend_label = NULL,
                        grid_res  = 300,
                        barwidth  = unit(10.5,   "cm"),
                        barheight = unit(0.75, "cm")) {
  ok  <- complete.cases(data.frame(lon, lat, z))
  fit <- MBA::mba.surf(data.frame(x = lon[ok], y = lat[ok], z = z[ok]),
                       no.X = grid_res, no.Y = grid_res, extend = TRUE)
  grd <- expand.grid(Longitude = fit$xyz.est$x, Latitude = fit$xyz.est$y)
  grd$z <- as.vector(fit$xyz.est$z)
  
  cb_guide <- guide_colorbar(title = legend_label, title.position = "top",
                             title.hjust = 0.5, barwidth = barwidth,
                             barheight = barheight, ticks.linewidth = 0.4)
  
  fscale <- if (pal == "RdBu") {
    scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
                         limits = limits, oob = scales::squish,
                         na.value = "transparent", guide = cb_guide)
  } else {
    scale_fill_viridis_c(option = pal, direction = -1, limits = limits,
                         oob = scales::squish, na.value = "transparent",
                         guide = cb_guide)
  }
  
  ggplot() +
    geom_raster(data = grd, aes(x = Longitude, y = Latitude, fill = z),
                interpolate = TRUE) +
    geom_sf(data = europe_sf, fill = NA, colour = "black", linewidth = 0.2) +
    fscale +
    coord_sf(xlim = range(lon[ok]), ylim = range(lat[ok]), expand = FALSE) +
    ggtitle(title) +
    theme_paper(leg = "right") +
    theme(legend.title = element_text(face = "bold", size = 18, hjust = 0.5),
          legend.text  = element_text(size = 9, face = "bold"),
          axis.text    = element_text(size = 9))
}

lon_u <- crdu[, 1]; lat_u <- crdu[, 2]

dyn_IS <- apply(Ysp_DYN_all[, , , 1L], c(1, 2), mean)
dyn_OS <- apply(Ysp_DYN_all[, , , 5L], c(1, 2), mean)
tr_IS  <- y_true_sp_all[, , 1L]
tr_OS  <- y_true_sp_all[, , 5L]

lims_IS <- lapply(seq_len(q), function(v)
  range(c(tr_IS[, v], dyn_IS[, v]), na.rm = TRUE))
lims_OS <- lapply(seq_len(q), function(v)
  range(c(tr_OS[, v], dyn_OS[, v]), na.rm = TRUE))

date_IS <- format(dates_vector[sp_is],      "%m-%Y")
date_OS <- format(dates_vector[sp_is + 24], "%m-%Y")

row_lbl <- function(txt)
  ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = txt,
           angle = 45, fontface = "italic", size = 6.5,
           hjust = 0.5, vjust = 0.5, colour = "grey20") +
  coord_cartesian(clip = "off") +
  theme_void()

make_col <- function(v, tr, dy, lims, lon, lat) {
  p_t <- geo_surface(lon, lat, tr[, v], title = NULL,
                     limits = lims[[v]], pal = var_palettes[v],
                     legend_label = var_labels[v])
  p_d <- geo_surface(lon, lat, dy[, v],
                     limits = lims[[v]], pal = var_palettes[v],
                     legend_label = var_labels[v])
  p_t / p_d   # flat: 2 righe, nient'altro
}

assemble_panel <- function(tr, dy, lims, tag, lon, lat) {
  cols   <- lapply(seq_len(q), make_col,
                   tr = tr, dy = dy, lims = lims, lon = lon, lat = lat)
  rl_col <- (row_lbl("True") / row_lbl("DYNBPS"))
  
  wrap_plots(c(list(rl_col), cols), nrow = 1,
             widths = c(0.25, rep(1, q))) +
    plot_layout(guides = "collect") +
    plot_annotation(title = tag,
                    theme = theme(plot.title = element_text(face  = "bold",
                                                            size  = 13,
                                                            hjust = 0.5,
                                                            margin = margin(b = 8)))) &
    theme(legend.position = "bottom")
}

fig7_IS <- assemble_panel(tr_IS, dyn_IS, lims_IS, "", lon_u, lat_u)
fig7_OS <- assemble_panel(tr_OS, dyn_OS, lims_OS, "", lon_u, lat_u)

ggsave("plots/fig7_IS_spatial.png", fig7_IS, width = 18, height = 9,  dpi = 300)
ggsave("plots/fig7_OS_spatial.png", fig7_OS, width = 18, height = 9,  dpi = 300)
cat("Fig 7 saved.\n")

# 11.5. interpolation figures ---------------------------------------------

snap_idx   <- round(seq(1, tmax, length.out = 4))
snap_dates <- format(dates_vector[snap_idx], "%m-%Y")

make_supp <- function(v) {
  z_tr   <- lapply(snap_idx, function(k) Ys[, v, k])
  z_dy   <- lapply(snap_idx, function(k) Y_postmean[, v, k])
  lims_v <- range(c(unlist(z_tr), unlist(z_dy)), na.rm = TRUE)
  n_snap <- length(snap_idx)
  
  # title = d: data numerica come intestazione di colonna sulla riga True
  r_true <- purrr::imap(snap_dates, function(d, k)
    geo_surface(lon_tr, lat_tr, z_tr[[k]], title = d,
                limits = lims_v, pal = var_palettes[v],
                legend_label = var_labels[v],
                barwidth = unit(20, "cm"), barheight = unit(0.6, "cm")))
  r_dyn  <- purrr::imap(snap_dates, function(d, k)
    geo_surface(lon_tr, lat_tr, z_dy[[k]],
                limits = lims_v, pal = var_palettes[v],
                legend_label = var_labels[v],
                barwidth = unit(20, "cm"), barheight = unit(0.6, "cm")))
  
  all_plots <- c(list(row_lbl("True")),   r_true,
                 list(row_lbl("DYNBPS")), r_dyn)
  
  wrap_plots(all_plots, nrow = 2, ncol = n_snap + 1L) +
    plot_layout(widths = c(0.20, rep(1, n_snap)), guides = "collect") &
    theme(legend.position = "bottom")
}

for (v in seq_len(q)) {
  fg <- make_supp(v)
  fn <- file.path("plots", paste0("figS", 11L + v, "_", tolower(var_labels[v])))
  ggsave(paste0(fn, ".png"), fg, width = 18, height = 8, dpi = 300)
  message("Saved S", 11L + v, " (", var_labels[v], ")")
}

cat("\nAll plots saved.\n")
