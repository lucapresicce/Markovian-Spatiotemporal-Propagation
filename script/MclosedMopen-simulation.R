###########################################################################
# M-CLOSED / M-OPEN SIMULATION + GRAPHICS 
###########################################################################

rm(list = ls()); gc()
Rcpp::sourceCpp("script/FFBS-DYNBPS-struct-v2.cpp")
source("script/FFBS-DYNBPS-v2.R")

library(mniw)
library(abind)
library(spFFBS)
library(spBPS)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)


# Global constants
tmax  <- 20L; tnew <- 5L
n     <- 100L; u   <- 50L
q     <- 3L;   p   <- 2L
BB    <- 50L
L     <- 200L
n_thr <- 5L

tau_true <- 0.80
phi_true <- 4.0
a_true   <- (1.0 / tau_true) - 1.0

Sigma <- matrix(c( 1.0, -0.3,  0.6,
                   -0.3,  1.2,  0.4,
                   0.6,  0.4,  1.0), q, q)


# 1. spatial structure & priors -------------------------------------------

set.seed(97)
coords_all <- matrix(runif((n + u) * 2L), ncol = 2L)
crd        <- coords_all[seq_len(n),     , drop = FALSE]
crdtilde   <- coords_all[n + seq_len(u), , drop = FALSE]
D          <- as.matrix(dist(crd))
D_u        <- as.matrix(dist(crdtilde))
D_us       <- spBPS:::arma_dist(rbind(crdtilde, crd))[seq_len(u), -seq_len(u)]
K_pr       <- exp(-1 * D)

m0       <- matrix(0, n + p, q)
C0       <- rbind(cbind(diag(0.005, p), matrix(0, p, n)),
                  cbind(matrix(0, n, p), K_pr))
prior_v2 <- make_prior_v2(m0, C0, nu0 = 3L, Psi0 = diag(q), p = p)

K_all  <- exp(-phi_true * as.matrix(dist(coords_all)))
W_big  <- rbind(cbind(diag(p),             matrix(0, p,     n + u)),
                cbind(matrix(0, n + u, p), K_all))
C0_big <- rbind(cbind(diag(0.005, p),      matrix(0, p,     n + u)),
                cbind(matrix(0, n + u, p), K_all))

set.seed(-100L)
theta0    <- mniw::rMNorm(1L, matrix(0, (n+u)+p, q), C0_big, Sigma)
theta_all <- array(0, c((n+u)+p, q, tmax+tnew))
for (tt in seq_len(tmax + tnew)) {
  noise            <- mniw::rMNorm(1L, matrix(0, (n+u)+p, q), W_big, Sigma)
  prev             <- if (tt == 1L) theta0 else theta_all[,, tt-1L]
  theta_all[,, tt] <- prev + noise
}


# 2. hyperparameters grids ------------------------------------------------

tau_C <- c(0.7, 0.8, 0.9)
phi_C <- c(2, 4, 6)
pg_C  <- spBPS:::expand_grid_cpp(rev(tau_C), rev(phi_C))
J_C   <- nrow(pg_C)
cat(sprintf("M-Closed grid: J=%d\n", J_C))

set.seed(2)
repeat {
  tau_O <- sort(runif(3L, 0.1, 1.0))
  phi_O <- sort(as.double(sample(setdiff(1:50, phi_true), 3L)))
  if (!any(abs(phi_O - phi_true) < 0.5)) break
}
pg_O <- spBPS:::expand_grid_cpp(rev(tau_O), rev(phi_O))
J_O  <- nrow(pg_O)
cat(sprintf("M-Open  grid: J=%d (fixed)\n  tau_O: %s\n  phi_O: %s\n",
            J_O,
            paste(round(tau_O, 3), collapse = ", "),
            paste(phi_O, collapse = ", ")))

extract_Sigma <- function(FF, wv, J) {
  indL <- sample.int(J, L, replace = TRUE, prob = wv)
  sp   <- vapply(seq_len(L), function(l)
    mniw::riwish(1L,
                 nu  = FF[[tmax]]$filtered_results[[indL[l]]]$nu,
                 Psi = FF[[tmax]]$filtered_results[[indL[l]]]$Psi),
    numeric(q * q))
  dim(sp) <- c(q, q, L); sp
}

threshold_w <- function(w, J) {
  wt <- pmax(w - 1.0/J, 0.0)
  s  <- sum(wt)
  if (s < 1e-12) return(rep(1.0/J, J))
  wt / s
}

proc_SI <- function(raw) {
  arr <- array(vapply(raw, as.numeric, numeric(2L*u*q)), c(2L*u, q, L))
  list(Y = arr[seq_len(u),,], Omega = arr[u + seq_len(u),,])
}


# 3. benchmarks computation -----------------------------------------------

compute_baseline_weights <- function(FF, J) {
  dens_list <- lapply(seq(2L, tmax), function(s) FF[[s]]$density_evaluations)
  
  # EW
  w_EW <- rep(1.0 / J, J)
  W_EW <- matrix(rep(w_EW, tmax - 1L), nrow = J)
  
  # Static local
  Wi_static <- matrix(0.0, n, J)
  log_dens <- do.call(rbind,
                      lapply(seq(2L, tmax), function(s) FF[[s]]$density_evaluations))
  
  w_raw    <- optimize_weights_proj_v2(exp(log_dens), lr = 1e-4, max_iter = 500L)
  w_Static_T <- w_raw / sum(w_raw)
  W_Static   <- matrix(rep(w_Static_T, tmax - 1L), nrow = J)

  # Top1
  W_Top1 <- matrix(0.0, J, tmax - 1L)
  for (s in seq_along(dens_list)) {
    j_stars     <- apply(dens_list[[s]], 1, which.max)
    W_Top1[, s] <- tabulate(j_stars, nbins = J) / n
    for(i in 1:J) { 
      W_Top1[i, s] <- ifelse(i == which.max(W_Top1[, s]), 1, 0) }
  }
  w_Top1_T <- W_Top1[, tmax - 1L]
  
  list(W_EW = W_EW, w_EW_T = w_EW,
       W_Static = W_Static, w_Static_T = w_Static_T,
       W_Top1 = W_Top1, w_Top1_T = w_Top1_T)
}

run_baseline <- function(FF, Wmat, wv_T, par_grid,
                         X_obs, Xtilde_arr, t_sp, t_out) {
  J <- length(wv_T)

  BS <- weighted_backward_sample_T_v2(
    G_beta = diag(p), rho = 1.0, D = D,
    ForwFilt = FF, L = as.integer(L),
    par_grid = par_grid, weights = wv_T,
    num_threads = n_thr
  )
  theta_post <- sapply(seq_len(tmax), function(tt) BS[[tt]], simplify = "array")
  Sigma_post <- extract_Sigma(FF, wv_T, J)
  
  TF     <- temporal_forecast_v2(
    G_beta = diag(p), rho = 1.0, D = D,
    par_grid = par_grid, ForwFilt = FF,
    X_all = X_obs, weights = Wmat,
    horiz = as.integer(tnew), L = as.integer(L),
    num_threads = n_thr
  )
  Y_pred <- abind::abind(TF, along = 4L)
  
  make_Wsp <- function(t_eff) {
    t_eff  <- min(t_eff, tmax)
    w_init <- matrix(1.0 / J, J, 1L)
    if (t_eff == 1L) return(w_init)
    cbind(w_init, Wmat[, seq_len(t_eff - 1L), drop = FALSE])
  }
  
  SI_in_raw <- spatial_interpolation_v2(
    G_beta = diag(p), rho = 1.0, Xu = Xtilde_arr[,, t_sp],
    D_s = D, D_u = D_u, D_us = D_us,
    par_grid = par_grid, ForwFilt = FF,
    weights = make_Wsp(t_sp),
    t = as.integer(t_sp), L = as.integer(L), num_threads = n_thr)
  SI_out_raw <- spatial_interpolation_v2(
    G_beta = diag(p), rho = 1.0, Xu = Xtilde_arr[,, t_out],
    D_s = D, D_u = D_u, D_us = D_us,
    par_grid = par_grid, ForwFilt = FF,
    weights = make_Wsp(min(t_out, tmax)),
    t = as.integer(t_out), L = as.integer(L), num_threads = n_thr)
  
  list(Sigma_post = Sigma_post, theta_post = theta_post,
       Y_pred = Y_pred,
       SI_in  = proc_SI(SI_in_raw),
       SI_out = proc_SI(SI_out_raw))
}


# 4. replications ---------------------------------------------------------

out_list <- vector("list", BB)
pb <- txtProgressBar(min = 0, max = BB, style = 3)

for (b in seq_len(BB)) {
  
  set.seed(b)
  cat(sprintf("\n=== Replication %d / %d ===\n", b, BB))
  
  # Data generation
  X_arr  <- array(runif((n+u)*p*(tmax+tnew)), c(n+u, p, tmax+tnew))
  Y_full <- array(0, c(n+u, q, tmax+tnew))
  for (tt in seq_len(tmax + tnew)) {
    Pt           <- cbind(X_arr[,,tt], diag(n+u))
    noise        <- mniw::rMNorm(1L, matrix(0, n+u, q), a_true * diag(n+u), Sigma)
    Y_full[,,tt] <- Pt %*% theta_all[,,tt] + noise
  }
  Y_obs      <- Y_full[seq_len(n),,   seq_len(tmax), drop = FALSE]
  X_obs      <- X_arr[ seq_len(n),,   ,              drop = FALSE]
  Xtilde_arr <- X_arr[ n+seq_len(u),, ,              drop = FALSE]
  t_sp  <- sample.int(tmax - 1L, 1L) + 1L
  t_out <- sample.int(tnew,      1L) + tmax
  cat("  Data generated.\n")
  
  # True model
  cat("  [ True      ] ... ")
  set.seed(b)
  out_true <- spFFBS_v2(
    Y = Y_obs, X = X_obs, D = D,
    grid  = list(tau = tau_true, phi = phi_true),
    prior = prior_v2, G_beta = NULL, rho = 1.0,
    do_BS = TRUE, do_forecast = TRUE, do_spatial = TRUE,
    L = L, tnew = tnew, X_all = X_obs,
    spatial = list(crd = crd, crdtilde = crdtilde,
                   Xtilde = Xtilde_arr[,, c(t_sp, t_out)],
                   t = c(t_sp, t_out)),
    n_threads = n_thr, verbose = FALSE
  )
  Sig_true <- vapply(seq_len(L), function(l)
    mniw::riwish(1L,
                 nu  = out_true$FF[[tmax]]$filtered_results[[1L]]$nu,
                 Psi = out_true$FF[[tmax]]$filtered_results[[1L]]$Psi),
    numeric(q*q)); dim(Sig_true) <- c(q, q, L)
  
  true_result <- list(
    Sigma_post = Sig_true,
    theta_post = sapply(seq_len(tmax),
                        function(tt) out_true$BS[[tt]], simplify = "array"),
    Y_pred = out_true$forecast$Y_pred,
    SI_in  = out_true$spatial[[1L]],
    SI_out = out_true$spatial[[2L]]
  )
  rm(out_true); gc()
  cat("done.\n")
  
  # M-CLOSED
  cat("  [ DYNBPS-C  ] ... ")
  set.seed(b)
  out_C <- spFFBS_v2(
    Y = Y_obs, X = X_obs, D = D,
    grid  = list(tau = tau_C, phi = phi_C),
    prior = prior_v2, G_beta = NULL, rho = 1.0,
    do_BS = TRUE, do_forecast = TRUE, do_spatial = TRUE,
    L = L, tnew = tnew, X_all = X_obs,
    spatial = list(crd = crd, crdtilde = crdtilde,
                   Xtilde = Xtilde_arr[,, c(t_sp, t_out)],
                   t = c(t_sp, t_out)),
    n_threads = n_thr, verbose = FALSE
  )
  
  Wg_thresh_C <- threshold_w(out_C$Wglobal, J_C)
  Sig_C <- extract_Sigma(out_C$FF, out_C$Wglobal, J_C)
  
  dynbps_C <- list(
    Sigma_post = Sig_C,
    theta_post = sapply(seq_len(tmax),
                        function(tt) out_C$BS[[tt]], simplify = "array"),
    Y_pred  = out_C$forecast$Y_pred,
    SI_in   = out_C$spatial[[1L]],
    SI_out  = out_C$spatial[[2L]],
    Wglobal = Wg_thresh_C,   # thresholded, for reporting
    Wi      = out_C$Wi        # raw, for diagnostics
  )
  cat("done.\n")
  
  cat("  [ Baselines-C] ... ")
  FF_C <- out_C$FF
  bw_C <- compute_baseline_weights(FF_C, J_C)
  rm(out_C); gc()
  
  EW_C     <- run_baseline(FF_C, bw_C$W_EW,    bw_C$w_EW_T,
                           pg_C, X_obs, Xtilde_arr, t_sp, t_out)
  Static_C <- run_baseline(FF_C, bw_C$W_Static, bw_C$w_Static_T,
                           pg_C, X_obs, Xtilde_arr, t_sp, t_out)
  Top1_C   <- run_baseline(FF_C, bw_C$W_Top1,   bw_C$w_Top1_T,
                           pg_C, X_obs, Xtilde_arr, t_sp, t_out)
  rm(FF_C); gc()
  cat("done.\n")
  
  closed_results <- list(
    DYNBPS = dynbps_C, EW = EW_C, Static = Static_C, Top1 = Top1_C, bw = bw_C)
  
  
  # M-OPEN
  cat("  [ DYNBPS-O  ] ... ")
  
  set.seed(b)
  tau_O <- sort(runif(3L, 0.1, 1.0))
  phi_O <- sort(round(runif(3L, 1, 50)))
  pg_O <- spBPS:::expand_grid_cpp(rev(tau_O), rev(phi_O))
  J_O  <- nrow(pg_O)
  cat(sprintf("M-Open  grid: J=%d (fixed)\n  tau_O: %s\n  phi_O: %s\n",
              J_O,
              paste(round(tau_O, 3), collapse = ", "),
              paste(phi_O, collapse = ", ")))
  
  set.seed(b)
  out_O <- spFFBS_v2(
    Y = Y_obs, X = X_obs, D = D,
    grid  = list(tau = tau_O, phi = phi_O),
    prior = prior_v2, G_beta = NULL, rho = 1.0,
    do_BS = TRUE, do_forecast = TRUE, do_spatial = TRUE,
    L = L, tnew = tnew, X_all = X_obs,
    spatial = list(crd = crd, crdtilde = crdtilde,
                   Xtilde = Xtilde_arr[,, c(t_sp, t_out)],
                   t = c(t_sp, t_out)),
    n_threads = n_thr, verbose = FALSE
  )
  
  Wg_thresh_O <- threshold_w(out_O$Wglobal, J_O)
  Sig_O <- extract_Sigma(out_O$FF, Wg_thresh_O, J_O)
  
  dynbps_O <- list(
    Sigma_post = Sig_O,
    theta_post = sapply(seq_len(tmax),
                        function(tt) out_O$BS[[tt]], simplify = "array"),
    Y_pred  = out_O$forecast$Y_pred,
    SI_in   = out_O$spatial[[1L]],
    SI_out  = out_O$spatial[[2L]],
    Wglobal = Wg_thresh_O,
    Wi      = out_O$Wi
  )
  cat("done.\n")
  
  cat("  [ Baselines-O] ... ")
  FF_O <- out_O$FF
  bw_O <- compute_baseline_weights(FF_O, J_O)
  rm(out_O); gc()
  
  EW_O     <- run_baseline(FF_O, bw_O$W_EW,    bw_O$w_EW_T,
                           pg_O, X_obs, Xtilde_arr, t_sp, t_out)
  Static_O <- run_baseline(FF_O, bw_O$W_Static, bw_O$w_Static_T,
                           pg_O, X_obs, Xtilde_arr, t_sp, t_out)
  Top1_O   <- run_baseline(FF_O, bw_O$W_Top1,   bw_O$w_Top1_T,
                           pg_O, X_obs, Xtilde_arr, t_sp, t_out)
  rm(FF_O); gc()
  cat("done.\n")
  
  open_results <- list(
    DYNBPS = dynbps_O, EW = EW_O, Static = Static_O, Top1 = Top1_O, bw = bw_O)
  
  out_list[[b]] <- list(
    True = true_result, Closed = closed_results, Open = open_results,
    Y_full = Y_full, t_sp = t_sp, t_out = t_out)
  
  setTxtProgressBar(pb, b)
}
close(pb)


# 5. graphics -------------------------------------------------------------

theta_true <- theta_all[seq_len(p),     , seq_len(tmax)]
omega_true <- theta_all[p + seq_len(n), , seq_len(tmax)]

# Palette and theme
SCH_COLORS <- c(
  "True"   = "#fdd835",
  "DYNBPS" = "#1976d2",
  "EW"     = "#c62828",
  "Static" = "#e65100",
  "Top1"   = "#2e7d32"
)
SCH_ORDER  <- names(SCH_COLORS)
SCH_LABELS <- c(
  "True"   = "True Model",
  "DYNBPS" = "DYNBPS",
  "EW"     = "Equal Weights",
  "Static" = "Static Weights",
  "Top1"   = "Top Model Selection"
)
SETTING_LABS <- c("Closed" = "M-Closed", "Open" = "M-Open")

theme_paper <- function(base = 18, leg = "none") {
  theme_bw(base_size = base) %+replace% theme(
    strip.text        = element_text(size = base, face = "bold"),
    strip.background  = element_rect(fill = "#eceff1", colour = "grey70"),
    axis.text         = element_text(size = base - 8, face = "bold"),
    axis.text.x       = element_text(size = base - 8, angle = 45, hjust = 0.5),
    axis.title        = element_text(size = base, face = "bold"),
    axis.title.x      = element_text(size = base - 2, face = "bold"),
    legend.text       = element_text(size = base - 2, face = "bold"),
    legend.title      = element_blank(),
    legend.key.size   = unit(1.3, "lines"),
    legend.position   = leg,
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "grey92")
  )
}

hline0   <- geom_hline(yintercept = 0,    linewidth = .4,
                       linetype = "dashed", colour = "grey40")
hline_95 <- geom_hline(yintercept = 0.95, linewidth = 1.1,
                       linetype = "dashed", colour = "firebrick")

scale_fill_sch  <- function()
  scale_fill_manual(values = SCH_COLORS, labels = SCH_LABELS)
scale_color_sch <- function()
  scale_colour_manual(values = SCH_COLORS, labels = SCH_LABELS)


# Data accessor
METHOD_TABLE <- bind_rows(
  data.frame(Setting="Closed", Scheme="True",   Setting_disp="M-Closed",
             stringsAsFactors=FALSE),
  data.frame(Setting="Open",   Scheme="True",   Setting_disp="M-Open",
             stringsAsFactors=FALSE),
  expand.grid(Setting=c("Closed","Open"),
              Scheme =c("DYNBPS","EW","Static","Top1"),
              stringsAsFactors=FALSE) %>%
    mutate(Setting_disp = SETTING_LABS[Setting])
)
get_res <- function(b, setting, scheme)
  if (scheme == "True") b$True else b[[setting]][[scheme]]


# 5.1. sigma posterior ----------------------------------------------------

cat("Building Sigma data...\n")

sigma_df <- local({
  rows <- list()
  for (ri in seq_len(nrow(METHOD_TABLE))) {
    sett  <- METHOD_TABLE$Setting[ri]
    sch   <- METHOD_TABLE$Scheme[ri]
    sdisp <- METHOD_TABLE$Setting_disp[ri]
    sp    <- sapply(out_list,
                    function(b) get_res(b,sett,sch)$Sigma_post,
                    simplify = "array")   # q×q×L×BB
    for (i in seq_len(q)) for (j in seq(i,q)) {
      samp  <- sp[i,j,,]
      pm    <- colMeans(samp)
      ci_lo <- apply(samp,2,quantile,.025)
      ci_hi <- apply(samp,2,quantile,.975)
      rows[[length(rows)+1]] <- data.frame(
        Setting=sdisp, Scheme=sch,
        Entry=paste0("(",i,",",j,")"),
        Bias    =mean(pm - Sigma[i,j]),
        Coverage=mean(Sigma[i,j] >= ci_lo & Sigma[i,j] <= ci_hi),
        SD      =mean(apply(samp,2,sd))
      )
    }
  }
  bind_rows(rows) %>%
    mutate(Scheme=factor(Scheme,SCH_ORDER),
           Entry=factor(Entry,
                        levels=c("(1,1)","(1,2)","(2,2)","(1,3)","(2,3)","(3,3)")))
})

sigma_xlabs <- c(
  expression(Sigma[1*","*1]), expression(Sigma[1*","*2]),
  expression(Sigma[2*","*2]), expression(Sigma[1*","*3]),
  expression(Sigma[2*","*3]), expression(Sigma[3*","*3])
)

frob_sigma_df <- local({
  rows <- list()
  for (ri in seq_len(nrow(METHOD_TABLE))) {
    sett  <- METHOD_TABLE$Setting[ri]
    sch   <- METHOD_TABLE$Scheme[ri]
    sdisp <- METHOD_TABLE$Setting_disp[ri]
    sp    <- sapply(out_list,
                    function(b) get_res(b,sett,sch)$Sigma_post,
                    simplify="array")
    norms <- vapply(seq_len(BB), function(r) {
      pm <- apply(sp[,,,r],c(1,2),mean)
      sqrt(sum((pm-Sigma)^2))
    }, numeric(1))
    rows[[length(rows)+1]] <- data.frame(Setting=sdisp,Scheme=sch,FrobNorm=norms)
  }
  bind_rows(rows) %>% mutate(Scheme=factor(Scheme,SCH_ORDER))
})

sigma_cov_new <- ggplot(
  sigma_df %>% mutate(Setting=factor(Setting,c("M-Closed","M-Open"))),
  aes(Entry, Coverage, fill=Scheme)) +
  geom_bar(stat="identity", position=position_dodge(.8), width=.75,
           key_glyph="rect",show.legend = F) +
  hline_95 +
  coord_cartesian(ylim=c(0,1.2)) +
  facet_wrap(~Setting, nrow=2) +
  scale_fill_manual(name=NULL, values=SCH_COLORS, labels=SCH_LABELS) +
  scale_x_discrete(labels=sigma_xlabs) +
  labs(x=NULL, y="Coverage (95%)") +
  theme_paper(base=20, leg="none")

sigma_frob_new <- ggplot(
  frob_sigma_df %>% mutate(Setting=factor(Setting,c("M-Closed","M-Open"))),
  aes(Scheme, FrobNorm, fill=Scheme)) +
  geom_violin(trim=TRUE, alpha=.85, key_glyph="rect") +
  geom_boxplot(width=.14, outlier.shape=NA, show.legend=FALSE) +
  coord_cartesian(ylim=c(0,1.5)) +
  facet_wrap(~Setting, nrow=2) +
  scale_fill_manual(name=NULL, values=SCH_COLORS, labels=SCH_LABELS) +
  scale_x_discrete(labels=SCH_LABELS[SCH_ORDER]) +
  labs(x=NULL, y=expression("||"*bar(Sigma)-Sigma[true]*"||"[F])) +
  theme_paper(base=20, leg="right")

sigma_new_plot <- (sigma_cov_new | sigma_frob_new) +
  plot_layout(guides="collect") &
  theme(legend.position="bottom")

png("plots/plot_sigma.png", width = 14, height = 8, units = "in", res = 300, family = "sans")
print(sigma_new_plot)
dev.off()


# 5.2. theta posterior -----------------------------------------------------

cat("Building Theta data...\n")

THETA_ENTRY_ORDER <- c("(1,1)","(1,2)","(1,3)","(2,1)","(2,2)","(2,3)")
theta_entry_labeller <- as_labeller(
  c("(1,1)"="Theta[1*','*1]","(1,2)"="Theta[1*','*2]","(1,3)"="Theta[1*','*3]",
    "(2,1)"="Theta[2*','*1]","(2,2)"="Theta[2*','*2]","(2,3)"="Theta[2*','*3]"),
  default=label_parsed)

theta_df <- local({
  rows <- list()
  for (ri in seq_len(nrow(METHOD_TABLE))) {
    sett  <- METHOD_TABLE$Setting[ri]
    sch   <- METHOD_TABLE$Scheme[ri]
    sdisp <- METHOD_TABLE$Setting_disp[ri]
    post  <- sapply(out_list,
                    function(b) get_res(b,sett,sch)$theta_post,
                    simplify="array")   # p×q×L×tmax×BB
    for (i in seq_len(p)) for (j in seq_len(q)) for (t in seq_len(tmax)) {
      samp  <- post[i,j,,t,]   # L×BB
      pm    <- colMeans(samp)
      ci_lo <- apply(samp,2,quantile,.025)
      ci_hi <- apply(samp,2,quantile,.975)
      tv    <- theta_true[i,j,t]
      rows[[length(rows)+1]] <- data.frame(
        Setting=sdisp, Scheme=sch,
        Entry=paste0("(",i,",",j,")"), Time=t,
        Bias    =mean(pm-tv),
        Coverage=mean(tv>=ci_lo & tv<=ci_hi),
        SD      =mean(apply(samp,2,sd)))
    }
  }
  bind_rows(rows) %>%
    mutate(Scheme=factor(Scheme,SCH_ORDER),
           Entry=factor(Entry,THETA_ENTRY_ORDER))
})

frob_theta_df <- local({
  rows <- list()
  for (ri in seq_len(nrow(METHOD_TABLE))) {
    sett  <- METHOD_TABLE$Setting[ri]
    sch   <- METHOD_TABLE$Scheme[ri]
    sdisp <- METHOD_TABLE$Setting_disp[ri]
    post  <- sapply(out_list,
                    function(b) get_res(b,sett,sch)$theta_post,
                    simplify="array")
    for (r in seq_len(BB)) for (t in seq_len(tmax)) {
      pm <- apply(post[seq_len(p),,,t,r],c(1,2),mean)
      rows[[length(rows)+1]] <- data.frame(
        Setting=sdisp,Scheme=sch,Replication=r,Time=t,
        FrobNorm=sqrt(sum((pm-theta_true[,,t])^2)))
    }
  }
  bind_rows(rows) %>% mutate(Scheme=factor(Scheme,SCH_ORDER))
})

med_theta <- frob_theta_df %>%
  group_by(Setting,Scheme,Time) %>%
  summarise(med=median(FrobNorm),.groups="drop")

make_theta_col <- function(sett_label) {
  
  cov_p <- ggplot(
    theta_df %>% filter(Setting == sett_label),
    aes(Time, Coverage, fill = Scheme, group = Scheme)) +
    geom_line(aes(colour = Scheme), linewidth = .9, key_glyph = "rect") +
    geom_point(aes(colour = Scheme), size = 1.2, show.legend = FALSE) +
    geom_hline(yintercept = .95, linewidth = .8,
               linetype = "dashed", colour = "firebrick") +
    facet_wrap(~ Entry, ncol = 3,
               labeller = labeller(Entry = theta_entry_labeller)) +
    scale_fill_manual(name = NULL, values = SCH_COLORS, labels = SCH_LABELS) +
    scale_colour_manual(guide = "none", values = SCH_COLORS) +
    scale_x_continuous(breaks = seq(1, tmax, by = 4)) +
    coord_cartesian(ylim = c(0, 1.25)) +
    labs(x = "", y = "Coverage (95%)",
         title = "") + #sett_label) +
    theme_paper(base = 22, leg = "bottom") +
    theme(plot.title = element_text(size = 36, face = "bold", hjust = .5))
  
  frob_p <- ggplot(
    frob_theta_df %>% filter(Setting == sett_label),
    aes(factor(Time), FrobNorm, fill = Scheme)) +
    geom_violin(trim = TRUE, scale = "width", alpha = .85,
                key_glyph = "rect") +
    geom_boxplot(width = .15, outlier.shape = NA,
                 position = position_dodge(.9),
                 colour = "grey30", show.legend = FALSE) +
    geom_point(
      data = med_theta %>% filter(Setting == sett_label),
      aes(factor(Time), med, colour = Scheme),
      size = 1.2, position = position_dodge(.5), show.legend = FALSE) +
    geom_line(
      data = med_theta %>% filter(Setting == sett_label),
      aes(factor(Time), med, colour = Scheme, group = Scheme),
      linewidth = .85, position = position_dodge(.5), show.legend = FALSE) +
    coord_cartesian(ylim = c(0, 2)) +
    scale_fill_manual(name = NULL, values = SCH_COLORS, labels = SCH_LABELS) +
    scale_colour_manual(guide = "none", values = SCH_COLORS) +
    scale_x_discrete(breaks = as.character(seq(1, tmax, by = 4))) +
    labs(x = "Time",
         y = expression("||" * bar(Theta) - Theta[true] * "||"[F])) +
    theme_paper(base = 22, leg = "bottom")
  
  cov_p / frob_p + plot_layout(heights = c(1, 1))
}

theta_col_C <- make_theta_col("M-Closed")
theta_col_O <- make_theta_col("M-Open")

png("plots/plot_theta_C.png", width = 14, height = 10, units = "in", res = 300, family = "sans")
print(theta_col_C)
dev.off()

png("plots/plot_theta_O.png", width = 14, height = 10, units = "in", res = 300, family = "sans")
print(theta_col_O)
dev.off()

# 5.3. omega posterior ----------------------------------------------------

cat("Building Omega data...\n")

omega_j_labeller <- as_labeller(
  setNames(paste0("Omega[",seq_len(q),"]"), as.character(seq_len(q))),
  default=label_parsed)

omega_df <- local({
  rows <- list()
  for (ri in seq_len(nrow(METHOD_TABLE))) {
    sett  <- METHOD_TABLE$Setting[ri]
    sch   <- METHOD_TABLE$Scheme[ri]
    sdisp <- METHOD_TABLE$Setting_disp[ri]
    # (p+n)×q×L×tmax×BB
    post  <- sapply(out_list,
                    function(b) get_res(b,sett,sch)$theta_post,
                    simplify="array")
    
    for (j in seq_len(q)) {
      for (t in seq_len(tmax)) {
        # n×L×BB
        samp <- post[p + seq_len(n), j, , t, , drop = FALSE]
        dim(samp) <- c(n, L, BB)
        
        tv <- omega_true[, j, t]   # n-vector
        
        # Vectorised: means and SDs across L (dim 2)
        means_nb <- apply(samp, c(1L, 3L), mean)   # n×BB
        sds_nb   <- apply(samp, c(1L, 3L), sd)     # n×BB
        
        # Coverage: quantiles across L — reshape to L×(n*BB), apply once
        samp_LnB <- aperm(samp, c(2L, 1L, 3L))               # L×n×BB
        samp_2d  <- matrix(samp_LnB, nrow = L, ncol = n*BB)  # L×(n*BB)
        q025     <- matrix(apply(samp_2d, 2L, quantile, .025), n, BB)
        q975     <- matrix(apply(samp_2d, 2L, quantile, .975), n, BB)
        tv_mat   <- matrix(tv, n, BB)
        cov_nb   <- (tv_mat >= q025) & (tv_mat <= q975)   # n×BB
        
        # Average over BB for each location → n rows
        rows[[length(rows)+1]] <- data.frame(
          Setting  = sdisp, Scheme = sch,
          i        = seq_len(n), j = j, Time = t,
          Bias     = rowMeans(means_nb) - tv,
          SD       = rowMeans(sds_nb),
          Coverage = rowMeans(cov_nb)
        )
      }
    }
  }
  bind_rows(rows) %>% mutate(Scheme=factor(Scheme,SCH_ORDER), j=factor(j))
})

frob_omega_df <- local({
  rows <- list()
  for (ri in seq_len(nrow(METHOD_TABLE))) {
    sett  <- METHOD_TABLE$Setting[ri]
    sch   <- METHOD_TABLE$Scheme[ri]
    sdisp <- METHOD_TABLE$Setting_disp[ri]
    post  <- sapply(out_list,
                    function(b) get_res(b,sett,sch)$theta_post,
                    simplify="array")
    for (r in seq_len(BB)) for (t in seq_len(tmax)) {
      pm <- apply(post[p+seq_len(n),,,t,r],c(1,2),mean)
      rows[[length(rows)+1]] <- data.frame(
        Setting=sdisp,Scheme=sch,Replication=r,Time=t,
        FrobNorm=sqrt(sum((pm-omega_true[,,t])^2)))
    }
  }
  bind_rows(rows) %>% mutate(Scheme=factor(Scheme,SCH_ORDER))
})

med_omega <- frob_omega_df %>%
  group_by(Setting,Scheme,Time) %>%
  summarise(med=median(FrobNorm),.groups="drop")

omega_line_df <- omega_df %>%
  group_by(Setting, Scheme, j, Time) %>%
  summarise(Coverage=mean(Coverage), .groups="drop") %>%
  mutate(Setting=factor(Setting,c("M-Closed","M-Open")),
         Scheme=factor(Scheme,SCH_ORDER))

make_omega_col <- function(sett_label) {
  
  cov_p <- ggplot(
    omega_line_df %>% filter(Setting == sett_label),
    aes(Time, Coverage, fill = Scheme, group = Scheme)) +
    geom_line(aes(colour = Scheme), linewidth = .9, key_glyph = "rect") +
    geom_point(aes(colour = Scheme), size = 1.2, show.legend = FALSE) +
    geom_hline(yintercept = .95, linewidth = .8,
               linetype = "dashed", colour = "firebrick") +
    facet_wrap(~ j, ncol = 1, labeller = omega_j_labeller) +
    scale_fill_manual(name = NULL, values = SCH_COLORS, labels = SCH_LABELS) +
    scale_colour_manual(guide = "none", values = SCH_COLORS) +
    scale_x_continuous(breaks = seq(1, tmax, by = 4)) +
    coord_cartesian(ylim = c(0, 1.25)) +
    labs(x = "", y = "Coverage (95%)") +
    theme_paper(base = 22, leg = "bottom")
  
  frob_p <- ggplot(
    frob_omega_df %>% filter(Setting == sett_label),
    aes(factor(Time), FrobNorm, fill = Scheme)) +
    geom_violin(trim = TRUE, scale = "width", alpha = .85,
                key_glyph = "rect") +
    geom_boxplot(width = .15, outlier.shape = NA,
                 position = position_dodge(.9),
                 colour = "grey30", show.legend = FALSE) +
    geom_point(
      data = med_omega %>% filter(Setting == sett_label),
      aes(factor(Time), med, colour = Scheme),
      size = 1.2, position = position_dodge(.5), show.legend = FALSE) +
    geom_line(
      data = med_omega %>% filter(Setting == sett_label),
      aes(factor(Time), med, colour = Scheme, group = Scheme),
      linewidth = .85, position = position_dodge(.5), show.legend = FALSE) +
    # coord_cartesian(ylim = c(6, 12)) +
    scale_fill_manual(name = NULL, values = SCH_COLORS, labels = SCH_LABELS) +
    scale_colour_manual(guide = "none", values = SCH_COLORS) +
    scale_x_discrete(breaks = as.character(seq(1, tmax, by = 4))) +
    labs(x = "Time",
         y = expression("||" * Omega - Omega[true] * "||"[F])) +
    theme_paper(base = 22, leg = "bottom")
  
  cov_p / frob_p + plot_layout(heights = c(1, 1))
}

omega_col_C <- make_omega_col("M-Closed")
omega_col_O <- make_omega_col("M-Open")

png("plots/plot_omega_C.png", width = 14, height = 10, units = "in", res = 300, family = "sans")
print(omega_col_C)
dev.off()

png("plots/plot_omega_O.png", width = 14, height = 10, units = "in", res = 300, family = "sans")
print(omega_col_O)
dev.off()


# 5.4. temporal forecast --------------------------------------------------

cat("Building Forecast data...\n")

y_fcst_true <- sapply(out_list,
                      function(b) b$Y_full[seq_len(n),,tmax+seq_len(tnew)],
                      simplify="array")   # n×q×tnew×BB

fcst_df <- local({
  rows <- list()
  for (ri in seq_len(nrow(METHOD_TABLE))) {
    sett  <- METHOD_TABLE$Setting[ri]
    sch   <- METHOD_TABLE$Scheme[ri]
    sdisp <- METHOD_TABLE$Setting_disp[ri]
    ypred <- sapply(out_list,
                    function(b) get_res(b,sett,sch)$Y_pred,
                    simplify="array")   # n×q×L×tnew×BB
    for (h in seq_len(tnew)) for (j in seq_len(q)) for (r in seq_len(BB)) {
      ps    <- ypred[,j,,tmax+h,r]
      tv    <- y_fcst_true[,j,h,r]
      pm    <- rowMeans(ps)
      ci_lo <- apply(ps,1,quantile,.025)
      ci_hi <- apply(ps,1,quantile,.975)
      rows[[length(rows)+1]] <- data.frame(
        Setting=sdisp, Scheme=sch, Horizon=h, Variable=j, Replication=r,
        RMSPE  =sqrt(mean((pm-tv)^2)),
        AbsBias=mean(abs(pm-tv)),
        PIWidth=mean(ci_hi-ci_lo))
    }
  }
  bind_rows(rows) %>%
    pivot_longer(c(RMSPE,AbsBias,PIWidth),
                 names_to="Metric",values_to="Value") %>%
    mutate(
      Scheme   = factor(Scheme,SCH_ORDER),
      Variable = factor(Variable, labels=paste0("Y[",seq_len(q),"]")),
      Metric   = factor(Metric,
                        levels=c("RMSPE","AbsBias","PIWidth"),
                        labels=c("RMSPE","Abs. Bias","PI Width")))
})

make_pred_row <- function(metric_label, y_label, x_label = "") {
  fcst_df %>%
    filter(Metric == metric_label) %>%
    mutate(Setting = factor(Setting, c("M-Closed","M-Open"))) %>%
    ggplot(aes(factor(Horizon), Value, fill = Scheme)) +
    geom_boxplot(width = .65, outlier.shape = NA,
                 position = position_dodge(.75)) +
    facet_grid(Setting ~ Variable, scales = "free_y",
               labeller = labeller(Variable = label_parsed)) +
    scale_fill_manual(values = SCH_COLORS, labels = SCH_LABELS, name = NULL) +
    labs(x = x_label, y = y_label) +
    theme_paper(base = 24, leg = "right")
}

pred_rmspe   <- make_pred_row("RMSPE",    "RMSPE",    x_label = "")
pred_piwidth <- make_pred_row("PI Width", "PI Width", x_label = "Forecast Horizon")

pred_new_plot <- (pred_rmspe / pred_piwidth) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

png("plots/plot_pred.png", width = 18, height = 12, units = "in", res = 300, family = "sans")
print(pred_new_plot)
dev.off()

