# -------------------------------------------------------------------
# AMORTIZED BAYESIAN FORECAST - SPATIOTEMPORAL TRANSFER LEARNING
# -------------------------------------------------------------------

# Requirements: keras (which uses tensorflow backend)
# Install if needed: install.packages("keras"); keras::install_keras()
library(keras)
library(tensorflow)

# Reproducibility
np <- import("numpy", convert = FALSE)
tf$random$set_seed(123L)
set.seed(123)


# Data generation ---------------------------------------------------------

# Problem dimensions & dataset
n <- 500          # locations
q <- 2            # variables per location
p <- 2            # number of covariates
t_steps <- 10     # time points
I <- 100          # number of independent samples (dataset size)
quantiles_k <- 3  # 50, 2.5, 97.5

# Data generator functions
# ALREADY DEFINED in ABI_genfun.R

# Build dataset of I couples
fixed <- tendency_gen(n = n, q = q, p = p, tmax = t_steps, phi = 4, seed = 42)
X_list <- vector("list", I)
Z_list <- vector("list", I)

# Define new dimensional parameters
qxp_X <- q + p + n   # input column dimension
qq_Z  <- q + q       # output column dimension

for (i in seq_len(I)) {
  Xi <- data_gen(n = n, q = q, p = p, tmax = t_steps, Theta = fixed$Theta)
  # assume Xi is now (n, qxp_X, t_steps)
  Zi <- calculate_post(Z = Xi, fixed = fixed, q = q, p = p)
  # assume Zi is now (n, qq_Z, quantiles_k)
  X_list[[i]] <- Xi
  Z_list[[i]] <- Zi
}

# Stack into arrays
X_all <- array(0, dim = c(I, n, qxp_X, t_steps))
Z_all <- array(0, dim = c(I, n, qq_Z, quantiles_k))

for (i in seq_len(I)) {
  X_all[i,,,] <- X_list[[i]][,,1:t_steps]
  Z_all[i,,,] <- Z_list[[i]]
}

cat("X_all shape:", paste(dim(X_all), collapse = " x "), "\n")  # I x n x (q+p+n) x t
cat("Z_all shape:", paste(dim(Z_all), collapse = " x "), "\n")  # I x n x (q+q) x 3

# Prepare inputs for Keras
# We'll reshape X -> (I, t, n*(q+p+n))
# And flatten targets Z -> (I, n*(q+q)*3)
I_train <- I  # use all for training

# Input tensor: (I, t, n*(q+p+n))
X_keras <- array(0, dim = c(I_train, t_steps, n * qxp_X))
for (i in 1:I_train) {
  tmp <- X_all[i,,,]  # (n, qxp_X, t)
  tmp2 <- aperm(tmp, c(3,1,2))  # (t, n, qxp_X)
  for (tt in 1:t_steps) {
    X_keras[i, tt, ] <- as.vector(tmp2[tt,,])
  }
}

# Target tensor: (I, n*(q+q)*3)
Z_keras <- array(0, dim = c(I_train, n * qq_Z * quantiles_k))
for (i in 1:I_train) {
  Z_keras[i, ] <- as.vector(Z_all[i,,,])
}

cat("X_keras shape:", paste(dim(X_keras), collapse = " x "), "\n")
cat("Z_keras shape:", paste(dim(Z_keras), collapse = " x "), "\n")

# Train/Val split
set.seed(97)
perm <- sample(1:I_train)
train_idx <- perm[1:round(0.8*I_train)]
val_idx   <- perm[(round(0.8*I_train)+1):I_train]
X_train <- X_keras[train_idx,,, drop = FALSE]
Y_train <- Z_keras[train_idx,, drop = FALSE]
X_val   <- X_keras[val_idx,,, drop = FALSE]
Y_val   <- Z_keras[val_idx,, drop = FALSE]
cat("Train size:", dim(X_train)[1], "Val size:", dim(X_val)[1], "\n")


# Amortized model ---------------------------------------------------------

# RNN: Encoder (GRU) -> latent -> predict Z quantiles
inputA <- layer_input(shape = c(t_steps, n * qxp_X), name = "inputA")

hA <- inputA %>%
  layer_gru(units = 64, return_sequences = FALSE, name = "gruA") %>%
  layer_dense(units = 128, activation = "relu", name = "denseA1")

# Output layer predicts flattened quantiles vector (n*(q+q)*3)
outA <- hA %>%
  layer_dense(units = n * qq_Z * quantiles_k, activation = "linear", name = "outA")

modelA <- keras_model(inputs = inputA, outputs = outA)
modelA$compile(optimizer = optimizer_adam(learning_rate = 1e-3), loss = "mse")

cat("Model summary:\n")
print(modelA$summary())

# Train RNN
cat("Training Recurrent Neural Network ...\n")
historyA <- modelA$fit(
  x = X_train, y = Y_train,
  validation_data = list(X_val, Y_val),
  epochs = 50L, steps_per_epoch = 10L,
  batch_size = 32L
)

historyA$history$loss |> plot(type = "l")
historyA$history$val_loss |> lines(col = 2)


# Predictive evaluation ---------------------------------------------------

# Generate new data
new_I <- 1

# Generate new input (like in your training loop)
X_new <- data_gen(n = n, tmax = t_steps, Theta = fixed$Theta)
# X_new has shape: (n, q+p+n, t_steps)

cat("New input X_new shape:", paste(dim(X_new), collapse = " x "), "\n")

# Prepare X_new for Keras
# We must reshape to match (batch_size, t_steps, n*(q+p+n))
# batch_size = new_I = 1
X_new_keras <- array(0, dim = c(new_I, t_steps, n * (q + p + n)))

tmp <- aperm(X_new, c(3,1,2))  # (t, n, q+p+n)
for (tt in 1:t_steps) {
  X_new_keras[1, tt, ] <- as.vector(tmp[tt,,])
}
cat("X_new_keras shape:", paste(dim(X_new_keras), collapse = " x "), "\n")

# Predict using the trained model
pred_new <- modelA$predict(X_new_keras)

# Reshape to (n, qq_Z, quantiles_k)
pred_new_reshaped <- array(pred_new, dim = c(n, (q+q), quantiles_k))
cat("Prediction shape:", paste(dim(pred_new_reshaped), collapse = " x "), "\n")

# Compare to ground truth
Z_true_new <- calculate_post(Z = X_new, fixed = fixed)

# Compute RMSE for the median quantile
pred_med_new <- pred_new_reshaped[,,1]
pred_low_new <- pred_new_reshaped[,,2]
pred_upp_new <- pred_new_reshaped[,,3]
true_med_new <- Z_true_new[,,1]

# ABF vs DYNBPS
rmse_new <- sqrt(mean((pred_med_new - true_med_new)^2))
cat(sprintf("RMSE median (new data): %.4f\n", rmse_new))

# ABF vs TRUE
rmse_abf <- sqrt(mean((pred_med_new[,1:2] - X_new[,1:2,11])^2))
cat(sprintf("RMSE median (new data): %.4f\n", rmse_abf))

# DYNBPS vs TRUE
rmse_bps <- sqrt(mean((true_med_new[,1:2] - X_new[,1:2,11])^2))
cat(sprintf("RMSE median (new data): %.4f\n", rmse_bps))

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
    labs(x = "Easting", y = "Northing", fill = bquote(hat(Y)[list(.(component),T+1)])) +
    coord_equal() +
    ggtitle(title) +
    theme_minimal()
}

# Y
plots_col1_Y <- list(
  plot_surface_interp(mat = X_new[,1:q,t_steps+1], component = 1, fixed_crd = fixed$crd, title = "True"),
  plot_surface_interp(mat = true_med_new, component = 1, fixed_crd = fixed$crd, title = "DYNBPS"),
  plot_surface_interp(mat = pred_med_new, component = 1, fixed_crd = fixed$crd, title = "Amortized")
)
plots_col2_Y <- list(
  plot_surface_interp(mat = X_new[,1:q,t_steps+1], component = 2, fixed_crd = fixed$crd, title = "True"),
  plot_surface_interp(mat = true_med_new, component = 2, fixed_crd = fixed$crd, title = "DYNBPS"),
  plot_surface_interp(mat = pred_med_new, component = 2, fixed_crd = fixed$crd, title = "Amortized")
)

row1_Y <- plot_grid(plotlist = plots_col1_Y, nrow = 1)
row2_Y <- plot_grid(plotlist = plots_col2_Y, nrow = 1)

# Save high-res PNG with Cairo
ggsave("plots/heatmap_amortized_Y.png", plot_grid(row1_Y, row2_Y, nrow = 2), dpi = 320, type = "cairo")


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
  
  if(component == 3 | component == 4) component <- component-2
  ggplot(df_grid, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdBu"))) + 
    labs(x = "Easting", y = "Northing", fill = bquote(hat(Omega)[list(.(component),T+1)])) +
    coord_equal() +
    ggtitle(title) +
    theme_minimal()
}

# Omega
plots_col1_Om <- list(
  plot_surface_interp(mat = fixed$Theta[-(1:p),,t_steps+1], component = 1, fixed_crd = fixed$crd, title = "True"),
  plot_surface_interp(mat = true_med_new, component = 3, fixed_crd = fixed$crd, title = "DYNBPS"),
  plot_surface_interp(mat = pred_med_new, component = 3, fixed_crd = fixed$crd, title = "Amortized")
  )
plots_col2_Om <- list(
  plot_surface_interp(mat = fixed$Theta[-(1:p),,t_steps+1], component = 2, fixed_crd = fixed$crd, title = "True"),
  plot_surface_interp(mat = true_med_new, component = 4, fixed_crd = fixed$crd, title = "DYNBPS"),
  plot_surface_interp(mat = pred_med_new, component = 4, fixed_crd = fixed$crd, title = "Amortized")
)

row1_Om <- plot_grid(plotlist = plots_col1_Om, nrow = 1)
row2_Om <- plot_grid(plotlist = plots_col2_Om, nrow = 1)

# Save high-res PNG with Cairo
ggsave("plots/heatmap_amortized_Om.png", plot_grid(row1_Om, row2_Om, nrow = 2), dpi = 320, type = "cairo")

