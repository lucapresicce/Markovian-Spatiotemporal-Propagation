#####################################################################
# AMORTIZED BAYESIAN FORECAST - DYNBPS vs MCMC TEACHER COMPARISON
#####################################################################

# Source order (assumed in working directory):
#   1. FFBS-DYNBPS-struct-v2.cpp   (Rcpp - compiled by ABF-genfun.R)
#   2. FFBS-DYNBPS-v2.R            (spFFBS_v2 wrapper)
#   3. ABF-genfun.R                (data gen + calculate_post + calculate_post_mcmc)
# -------------------------------------------------------------------

library(keras)
library(tensorflow)
library(cowplot)
library(ggplot2)
library(patchwork)
library(magick)

source("script_finals/ABF-genfun.R")

# Reproducibility
py_require("tensorflow")
np <- import("numpy", convert = FALSE)
tf$random$set_seed(123L)
set.seed(123)


# 0. parameters -----------------------------------------------------------

n        <- 250L   # spatial locations
q        <- 2L     # variables per location
p        <- 2L     # covariates (intercept + 1)
t_steps  <- 10L    # time points for fitting; T+1 is prediction target
I        <- 50L    # training instances
n_test   <- 10L    # test instances for evaluation
L        <- 200L   # posterior samples per instance
N_THREADS <- 4L    # threads for spFFBS_v2 (DYNBPS only)

# MCMC settings
N_ITER_MC   <- 2000L
N_BURNIN_MC <- 1000L

# Derived Keras dimensions
qxp_X      <- q + p + n   # input feature dimension
qq_Z       <- q + q        # output, each q columns
quantiles_k <- 3L           # 50th / 2.5th / 97.5th


# 1. fixed spatial structure ----------------------------------------------

fixed <- tendency_gen(n = n, q = q, p = p, tmax = t_steps, phi = 4, seed = 42)
cat(sprintf("Fixed structure generated.  n=%d  q=%d  p=%d  T=%d\n",
            n, q, p, t_steps))


# 2. generate datasets ----------------------------------------------------

cat(sprintf("\nGenerating %d datasets...\n", I))
set.seed(42)
X_list <- vector("list", I)
for (i in seq_len(I)) {
  X_list[[i]] <- data_gen(n = n, q = q, p = p,
                           tmax = t_steps, Theta = fixed$Theta)
}
cat("  Done.\n")


# 3. DYNBPS posteriors ----------------------------------------------------

cat(sprintf("\n[Teacher 1 / DYNBPS]  Computing posteriors for %d datasets...\n", I))
time_bps_each <- numeric(I)
Z_list_bps    <- vector("list", I)

for (i in seq_len(I)) {
  t0 <- proc.time()
  Z_list_bps[[i]] <- calculate_post(X_list[[i]], fixed, q = q, p = p,
                                     n_threads = N_THREADS)
  time_bps_each[i] <- (proc.time() - t0)["elapsed"]
  if (i %% 10L == 0L)
    cat(sprintf("  [BPS] %d/%d done  (last: %.1f sec)\n", i, I, time_bps_each[i]))
}

time_bps_mean <- mean(time_bps_each)
cat(sprintf("  DYNBPS: %.2f sec/instance  (total: %.1f min)\n",
            time_bps_mean, sum(time_bps_each) / 60))


# 4. MCMC posteriors ------------------------------------------------------

cat(sprintf("\n[Teacher 2 / MCMC]  Computing posteriors for %d datasets...\n", I))
cat(sprintf("  Settings: n_iter=%d  n_burnin=%d  L=%d\n",
            N_ITER_MC, N_BURNIN_MC, L))
cat("  Expected time: ~", round(I * 180 / 60), "min (estimate based on n=100, T=10)\n")

time_mc_each <- numeric(I)
Z_list_mc    <- vector("list", I)

for (i in seq_len(I)) {
  t0 <- proc.time()
  Z_list_mc[[i]] <- calculate_post_mcmc(
    Z              = X_list[[i]],
    fixed          = fixed,
    q              = q,  p = p,
    n_iter         = N_ITER_MC,
    n_burnin       = N_BURNIN_MC,
    L              = L,
    phi_prior_shape = 2,
    phi_prior_rate  = 0.5
  )
  time_mc_each[i] <- (proc.time() - t0)["elapsed"]
  if (i %% 1L == 0L)
    cat(sprintf("  [MC]  %d/%d done  (last: %.1f sec  acc: %.2f)\n",
                i, I, time_mc_each[i],
                attr(Z_list_mc[[i]], "acc_rate")))
}

time_mc_mean <- mean(time_mc_each)
cat(sprintf("  MCMC: %.2f sec/instance  (total: %.1f min)\n",
            time_mc_mean, sum(time_mc_each) / 60))


# 4b. saving teacher-supervised data --------------------------------------

supervised_data <- list(
  X_list     = X_list,
  Z_list_bps = Z_list_bps,
  Z_list_mc  = Z_list_mc,
  settings = list(
    I = I,
    n = n,
    q = q,
    p = p,
    t_steps = t_steps,
    fixed = fixed,
    N_ITER_MC = N_ITER_MC,
    N_BURNIN_MC = N_BURNIN_MC,
    L = L
  ),
  timing = list(
    time_bps_each = time_bps_each,
    time_bps_mean = time_bps_mean,
    time_mc_each  = time_mc_each,
    time_mc_mean  = time_mc_mean
  )
)

saveRDS(
  supervised_data,
  file = "data/supervised_data.rds",
  compress = "xz"
)

# results <- readRDS("data/supervised_data.rds")
# X_list     <- results$X_list
# Z_list_bps <- results$Z_list_bps
# Z_list_mc  <- results$Z_list_mc

# 5. timing summary -------------------------------------------------------

cat("\n=== Teacher computational cost ===\n")
cat(sprintf("  DYNBPS  : %.2f sec/instance\n", time_bps_mean))
cat(sprintf("  MCMC    : %.2f sec/instance\n", time_mc_mean))
cat(sprintf("  Ratio   : %.1fx  (MCMC / DYNBPS)\n", time_mc_mean / time_bps_mean))


# 6. build keras arrays ---------------------------------------------------

make_X_keras <- function(X_list_, I_, n_, t_, qxp_) {
  out <- array(0, dim = c(I_, t_, n_ * qxp_))
  for (i in seq_len(I_)) {
    tmp <- aperm(X_list_[[i]][,, 1:t_, drop = FALSE], c(3, 1, 2))  # t x n x qxp
    for (tt in seq_len(t_)) out[i, tt, ] <- as.vector(tmp[tt,,])
  }
  out
}

make_Z_keras <- function(Z_list_, I_, n_, qq_, k_) {
  out <- array(0, dim = c(I_, n_ * qq_ * k_))
  for (i in seq_len(I_)) out[i,] <- as.vector(Z_list_[[i]])
  out
}

X_keras     <- make_X_keras(X_list, I, n, t_steps, qxp_X)
Z_keras_bps <- make_Z_keras(Z_list_bps, I, n, qq_Z, quantiles_k)
Z_keras_mc  <- make_Z_keras(Z_list_mc,  I, n, qq_Z, quantiles_k)

cat("\nArray shapes:\n")
cat("  X_keras     :", paste(dim(X_keras),     collapse = " x "), "\n")
cat("  Z_keras_bps :", paste(dim(Z_keras_bps), collapse = " x "), "\n")
cat("  Z_keras_mc  :", paste(dim(Z_keras_mc),  collapse = " x "), "\n")


# 7. train / val split ----------------------------------------------------

set.seed(97)
perm      <- sample(seq_len(I))
train_idx <- perm[1:round(0.8 * I)]
val_idx   <- perm[(round(0.8 * I) + 1):I]

X_train_k <- X_keras[train_idx,,, drop = FALSE]
X_val_k   <- X_keras[val_idx,,,  drop = FALSE]

Y_train_bps <- Z_keras_bps[train_idx,, drop = FALSE]
Y_val_bps   <- Z_keras_bps[val_idx,,  drop = FALSE]

Y_train_mc  <- Z_keras_mc[train_idx,, drop = FALSE]
Y_val_mc    <- Z_keras_mc[val_idx,,   drop = FALSE]

cat(sprintf("\nTrain: %d  Val: %d\n", length(train_idx), length(val_idx)))


# 8. GRU architecture -----------------------------------------------------

build_gru_model <- function(t_, n_, qxp_, qq_, k_) {
  inp <- layer_input(shape = c(t_, n_ * qxp_))
  h   <- inp %>%
    layer_gru(units = 64L, return_sequences = FALSE) %>%
    layer_dense(units = 128L, activation = "relu")
  out <- h %>%
    layer_dense(units = n_ * qq_ * k_, activation = "linear")
  mdl <- keras_model(inputs = inp, outputs = out)
  mdl$compile(optimizer = optimizer_adam(learning_rate = 1e-3), loss = "mse")
  mdl
}


# 9. train model A  (DYNBPS teacher) --------------------------------------

cat("\n--- Training GRU - DYNBPS teacher ---\n")
set.seed(123); tf$random$set_seed(123L)
modelA <- build_gru_model(t_steps, n, qxp_X, qq_Z, quantiles_k)
cat("Model parameters:", format(modelA$count_params(), big.mark = ","), "\n")

historyA <- modelA$fit(
  x                = X_train_k,
  y                = Y_train_bps,
  validation_data  = list(X_val_k, Y_val_bps),
  epochs           = 100L,
  steps_per_epoch  = 32L,
  batch_size       = 32L,
  verbose          = 1L
)


# 10. train model B  (MCMC teacher) ---------------------------------------

cat("\n--- Training GRU - MCMC teacher ---\n")
set.seed(123); tf$random$set_seed(123L)
modelB <- build_gru_model(t_steps, n, qxp_X, qq_Z, quantiles_k)

historyB <- modelB$fit(
  x                = X_train_k,
  y                = Y_train_mc,
  validation_data  = list(X_val_k, Y_val_mc),
  epochs           = 100L,
  steps_per_epoch  = 32L,
  batch_size       = 32L,
  verbose          = 1L
)


# 11. evaluation on test datasets -----------------------------------------

cat(sprintf("\n--- Generating %d test datasets ---\n", n_test))
set.seed(1995)
X_test_list <- vector("list", n_test)
Z_test_bps  <- vector("list", n_test)
Z_test_mc   <- vector("list", n_test)

for (i in seq_len(n_test)) {
  X_test_list[[i]] <- data_gen(n = n, q = q, p = p,
                               tmax = t_steps, Theta = fixed$Theta)
  Z_test_bps[[i]]  <- calculate_post(X_test_list[[i]], fixed,
                                     q = q, p = p, n_threads = N_THREADS)
  Z_test_mc[[i]]   <- calculate_post_mcmc(X_test_list[[i]], fixed,
                                          q = q, p = p,
                                          n_iter   = N_ITER_MC,
                                          n_burnin = N_BURNIN_MC,
                                          L        = L)
  cat(sprintf("  test %d/%d done\n", i, n_test))
}

# Ground truth
Y_true_all  <- sapply(X_test_list, function(Xi) Xi[, 1:q, t_steps + 1L],
                      simplify = "array")
Om_true_all <- replicate(n_test, fixed$Theta[(p+1L):(p+n),, t_steps + 1L])

# Metric helpers
rmspe_q    <- function(c, y) {sqrt(mean((c - y)^2))}
coverage_q <- function(lo, hi, y) {mean(y >= lo & y <= hi)}
IS95_q     <- function(lo, hi, y, a = 0.05) {
  mean((hi - lo) + (2/a) * (pmax(lo - y, 0) + pmax(y - hi, 0)))
}

results <- data.frame()
for (i in seq_len(n_test)) {
  Xi      <- X_test_list[[i]]
  Y_true  <- Y_true_all[,, i]
  Om_true <- Om_true_all[,, i]
  
  X_k <- array(0, dim = c(1L, t_steps, n * qxp_X))
  tmp <- aperm(Xi[,, 1:t_steps, drop = FALSE], c(3, 1, 2))
  for (tt in seq_len(t_steps)) X_k[1, tt,] <- as.vector(tmp[tt,,])
  
  pred_A <- array(modelA$predict(X_k, verbose = 0L), c(n, qq_Z, quantiles_k))
  pred_B <- array(modelB$predict(X_k, verbose = 0L), c(n, qq_Z, quantiles_k))
  
  for (tgt in c("Y", "Omega")) {
    cols <- if (tgt == "Y") 1:q else (q + 1):(2 * q)
    y_t  <- if (tgt == "Y") Y_true else Om_true
    for (mth in c("DYNBPS", "MCMC", "NN_DYNBPS", "NN_MCMC")) {
      W_m <- switch(mth,
                    DYNBPS   = Z_test_bps[[i]],
                    MCMC     = Z_test_mc[[i]],
                    NN_DYNBPS = pred_A,
                    NN_MCMC  = pred_B)
      results <- rbind(results, data.frame(
        Dataset  = i, Method = mth, Target = tgt,
        RMSPE    = rmspe_q(W_m[, cols, 1], y_t),
        Coverage = coverage_q(W_m[, cols, 2], W_m[, cols, 3], y_t),
        IS95     = IS95_q(W_m[, cols, 2], W_m[, cols, 3], y_t)
      ))
    }
  }
}


# 12. metrics table -------------------------------------------------------

agg <- aggregate(cbind(RMSPE, Coverage, IS95) ~ Method + Target,
                 data = results, FUN = mean)

agg_Y  <- agg[agg$Target == "Y",     c("Method", "RMSPE", "Coverage", "IS95")]
agg_Om <- agg[agg$Target == "Omega", c("Method", "RMSPE", "Coverage", "IS95")]
names(agg_Y)  <- c("Method", "Y_RMSPE",  "Y_Cov95",  "Y_IS95")
names(agg_Om) <- c("Method", "Om_RMSPE", "Om_Cov95", "Om_IS95")

tbl <- merge(agg_Y, agg_Om, by = "Method")

# Teacher cost per training instance
tbl$Time_sec <- c(DYNBPS    = time_bps_mean,
                  MCMC      = time_mc_mean,
                  NN_DYNBPS = NA_real_,
                  NN_MCMC   = NA_real_)[tbl$Method]

row_order <- c("DYNBPS", "MCMC", "NN_DYNBPS", "NN_MCMC")
tbl <- tbl[match(row_order, tbl$Method),]
tbl[, -1] <- round(tbl[, -1], 4)

cat("\n")
cat(sprintf("%-12s | %6s %7s %6s | %8s %7s %6s | %9s\n",
            "Method", "RMSPE", "Cov95", "IS95",
            "Om_RMSPE", "Cov95", "IS95", "Time(sec)"))
cat(strrep("-", 72), "\n")
for (r in seq_len(nrow(tbl))) {
  cat(sprintf("%-12s | %6.4f %7.4f %6.4f | %8.4f %7.4f %6.4f | %9s\n",
              tbl$Method[r],
              tbl$Y_RMSPE[r], tbl$Y_Cov95[r], tbl$Y_IS95[r],
              tbl$Om_RMSPE[r], tbl$Om_Cov95[r], tbl$Om_IS95[r],
              ifelse(is.na(tbl$Time_sec[r]), "-",
                     sprintf("%.1f", tbl$Time_sec[r]))))
}
cat(strrep("-", 72), "\n")

cat("\n=== Teacher cost per training instance ===\n")
cat(sprintf("  DYNBPS : %.2f sec\n", time_bps_mean))
cat(sprintf("  MCMC   : %.2f sec  (%.1fx slower)\n",
            time_mc_mean, time_mc_mean / time_bps_mean))


# 13. surface plots -------------------------------------------------------

# Use last test dataset for illustration
i_plt  <- n_test
Xi_plt <- X_test_list[[i_plt]]
X_k_plt <- array(0, c(1L, t_steps, n * qxp_X))
tmp_plt <- aperm(Xi_plt[,, 1:t_steps, drop = FALSE], c(3,1,2))
for (tt in seq_len(t_steps)) X_k_plt[1, tt,] <- as.vector(tmp_plt[tt,,])

pred_A_plt <- array(modelA$predict(X_k_plt, verbose=0L), c(n, qq_Z, quantiles_k))
pred_B_plt <- array(modelB$predict(X_k_plt, verbose=0L), c(n, qq_Z, quantiles_k))
W_bps_plt  <- Z_test_bps[[i_plt]]

# plot_surface_mba is defined in ABF-genfun.R
pmba <- function(mat, comp, title, limits, leg = NULL)
  plot_surface_mba(mat, fixed$crd,
                   component    = comp,
                   title        = title,
                   limits       = limits,
                   legend_label = leg)

cth <- theme(legend.position = "right",
             legend.title    = element_text(),
             legend.text     = element_text(size = 16, face = "bold"))

lims_y1 <- range(Xi_plt[, 1, t_steps + 1L],
                  pred_A_plt[, 1, 1], pred_B_plt[, 1, 1], na.rm = TRUE)
lims_y2 <- range(Xi_plt[, 2, t_steps + 1L],
                  pred_A_plt[, 2, 1], pred_B_plt[, 2, 1], na.rm = TRUE)

lab_Y1  <- expression(hat(Y)["1, T+1"])
lab_Y2  <- expression(hat(Y)["2, T+1"])

# Plots Y
plots_Y <- list(
  pmba(Xi_plt[, 1:q, t_steps + 1L], 1, "True",        lims_y1, lab_Y1),
  pmba(Xi_plt[, 1:q, t_steps + 1L], 2, "True",        lims_y2, lab_Y2),
  pmba(pred_A_plt[, 1:q, 1],        1, "Amortized (DYNBPS)",  lims_y1, lab_Y1),
  pmba(pred_A_plt[, 1:q, 1],        2, "Amortized (DYNBPS)",  lims_y2, lab_Y2),
  pmba(pred_B_plt[, 1:q, 1],        1, "Amortized (MCMC)",    lims_y1, lab_Y1),
  pmba(pred_B_plt[, 1:q, 1],        2, "Amortized (MCMC)",    lims_y2, lab_Y2)
)

row1   <- (plots_Y[[1]] | plots_Y[[3]] | plots_Y[[5]]) +
  plot_layout(guides = "collect") & cth
row2   <- (plots_Y[[2]] | plots_Y[[4]] | plots_Y[[6]]) +
  plot_layout(guides = "collect") & cth
grid_Y <- row1 / row2

lims_Om1 <- range(fixed$Theta[(p+1L):(p+n), 1, t_steps+1L],
                  pred_A_plt[, 3, 1], pred_B_plt[, 3, 1], na.rm = TRUE)
lims_Om2 <- range(fixed$Theta[(p+1L):(p+n), 2, t_steps+1L],
                  pred_A_plt[, 4, 1], pred_B_plt[, 4, 1], na.rm = TRUE)

lab_Om1 <- expression(hat(Omega)["1, T+1"])
lab_Om2 <- expression(hat(Omega)["2, T+1"])

Om_true <- fixed$Theta[(p+1L):(p+n), , t_steps+1L]
Om_A    <- pred_A_plt[, (q+1):(2*q), 1]
Om_B    <- pred_B_plt[, (q+1):(2*q), 1]

# Plots Omega
plots_Om <- list(
  pmba(Om_true, 1, "True",        lims_Om1, lab_Om1),
  pmba(Om_true, 2, "True",        lims_Om2, lab_Om2),
  pmba(Om_A,    1, "Amortized (DYNBPS)",  lims_Om1, lab_Om1),
  pmba(Om_A,    2, "Amortized (DYNBPS)",  lims_Om2, lab_Om2),
  pmba(Om_B,    1, "Amortized (MCMC)",    lims_Om1, lab_Om1),
  pmba(Om_B,    2, "Amortized (MCMC)",    lims_Om2, lab_Om2)
)

row1    <- (plots_Om[[1]] | plots_Om[[3]] | plots_Om[[5]]) +
  plot_layout(guides = "collect") & cth
row2    <- (plots_Om[[2]] | plots_Om[[4]] | plots_Om[[6]]) +
  plot_layout(guides = "collect") & cth
grid_Om <- row1 / row2

ggsave("plots/heatmap_amortized_Y.png",    grid_Y,  width=18, height=8,  dpi=300)
ggsave("plots/heatmap_amortized_Om.png",   grid_Om, width=18, height=8,  dpi=300)

# image cut
trim_png <- function(path, border = "20x20") {
  image_read(path) |>
    image_trim() |>
    image_border(color = "white", geometry = border) |>
    image_write(path)
}
trim_png("plots/heatmap_amortized_Y.png")
trim_png("plots/heatmap_amortized_Om.png")


