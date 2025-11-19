## WEIGHTS DYNAMIC - SIMULATION ###########################################

# Working packages
library(spBPS)
library(mniw)
library(classInt)
library(RColorBrewer)
library(sp)
library(fields)
library(parallel)
library(foreach)
library(doParallel)

library(dplyr)
library(tidyr)
library(ggplot2)
library(deldir)
library(scales)
library(patchwork)
# 
# # Working cpp functions
# Rcpp::sourceCpp("Rcpp/FFBS.cpp")
# Rcpp::sourceCpp("Rcpp/FFBS-BPS.cpp")
# Rcpp::sourceCpp("Rcpp/FFBS-parallel.cpp")
# Rcpp::sourceCpp("Rcpp/FFBS-prediction.cpp")
# Rcpp::sourceCpp("Rcpp/FFBS-parallel3.cpp")
# Rcpp::sourceCpp("Rcpp/FFBS-parallel-weights.cpp")
# 
# # Working R function
# compute_Wt <- function(out_FF, t, n, tau_seq, phi_seq) {
#   
#   cat("Computing Dynamic Bayesian Predictive Stacking Weights ...\n")
#   tictoc::tic()
#   
#   # Parallel computing over n
#   Wi <- foreach(i = 1:n, .combine = "rbind", .packages = "spBPS") %dopar% {
#     epd_i <- sapply(2:t, function(a) out_FF[[a]]$density_evaluations[i, ])
#     Wi_i <- spBPS::conv_opt(scores = exp(t(epd_i)))
#     # Wi_i <- optimize_weights_proj(scores = exp(t(epd_i)))
#     t(Wi_i)  # Ensure row structure
#   }
#   J <- ncol(Wi)
#   colnames(Wi) <- paste0(1:J, " - ", rev(tau_seq), ", " , rev(phi_seq))
#   tictoc::toc()
#   
#   # Weights
#   return(Wi)
# }

# Associated R package
# devtools::install_github("lucapresicce/spFFBS") # if not previously installed
library(spFFBS)
library(spBPS)

# Data Simulation ---------------------------------------------------------

# Synthetic data dimensions
tmax <- 100
n    <- 500
q    <- 3
p    <- 2

# Coordinates and true parameters
crd <- matrix(runif((n)*2), (n), 2)
D   <- as.matrix(dist(crd))
phi <- 4
K   <- exp(-phi*D)
W <- rbind( cbind(diag(p), matrix(0, p, (n))), cbind(matrix(0, (n), p), K) )
tau <- 0.80
a <- ((1/tau)-1); V <- a*diag(n)
Sigma  <- matrix(c(1, -0.3, 0.6, -0.3, 1.2, 0.4, 0.6, 0.4, 1), q, q)

# Prior information
m0     <- matrix(0, n+p, q)
C0 <- rbind( cbind(diag(0.005, p), matrix(0, p, (n))), cbind(matrix(0, (n), p), K) )

# Initial state
theta0 <- mniw::rMNorm(n = 1, Lambda = m0, SigmaR = C0, SigmaC = Sigma)

# Predictors and containers
X <- array(0, c(n, p, tmax))
P <- array(0, c(n, n+p, tmax))
G <- array(0, c(n+p, n+p, tmax))
Y <- array(0, c(n, q, tmax))
theta <- array(0, c(n+p, q, tmax))

# Generate data from the model
set.seed(4-8-15-16-23-42)
for (t in 1:(tmax)) {
  if (t >= 2) {  
    
    X[,,t]     <- matrix(runif(n*p), n, p)
    P[,,t]     <- cbind(X[,,t], diag(n))
    G[,,t]     <- diag(p+n)
    theta[,,t] <- G[,,t] %*% theta[,,t-1] + mniw::rMNorm(n = 1, Lambda = m0, SigmaR = W, SigmaC = Sigma)
    Y[,,t]     <- P[,,t] %*% theta[,,t] + mniw::rMNorm(n = 1, Lambda = matrix(0, n, q), SigmaR = V, SigmaC = Sigma)
  } else {
    
    X[,,t]     <- matrix(runif(n*p), n, p)
    P[,,t]     <- cbind(X[,,t], diag(n))
    G[,,t]     <- diag(p+n)
    theta[,,t] <- G[,,t] %*% theta0 + mniw::rMNorm(n = 1, Lambda = m0, SigmaR = W, SigmaC = Sigma)
    Y[,,t]     <- P[,,t] %*% theta[,,t] + mniw::rMNorm(n = 1, Lambda = matrix(0, n, q), SigmaR = V, SigmaC = Sigma)
  }
}

# Prior information
m0     <- matrix(0, n+p, q)
C0 <- rbind( cbind(diag(0.005, p), matrix(0, p, n)), cbind(matrix(0, n, p), K) )
nu0 <- 3
Psi0 <- diag(q)

# Setting Bayesian Predictive Stacking ------------------------------------

# Define competetive (J) models
tau_seq <- c(0.7, 0.8, 0.9)
phi_seq <- c(2, 4, 6)
par_grid <- spBPS::expand_grid_cpp(rev(tau_seq), rev(phi_seq))
J <- nrow(par_grid)

# prior information
prior <- list("m" = m0, "C" = C0, "nu" = nu0, "Psi" = Psi0)


# DYNBPS forward filtering backward sampling ------------------------------

out <- spFFBS::spFFBS(Y = Y, G = G, P = P, D = D,
                      grid = list(tau = tau_seq, phi = phi_seq), prior = prior, L = 200)

# DYNBPS weights computation ----------------------------------------------

# Compute weights for all time points
Wtemp <- array(1/J, c(n, J, tmax))
tictoc::tic()
for (t in 2:tmax) {
  density_list <- lapply(2:t, function(tt) out$FF[[tt]]$density_evaluations)
  Wi <- spFFBS:::compute_Wt_cpp(density_list, n = n, t = t, lr = 0.05, max_iter = 500, n_threads = 1)
  Wtemp[,,t] <- Wi
}
tictoc::toc()

# global weights
(Wglobal <- colMeans(Wi) |> matrix())

# consensus weights
pref_model <- apply(Wi, 1, which.max)
Wcons <- matrix(table(pref_model)/sum(table(pref_model)))
if(length(Wcons)!=J) { 
  Wcons0 <-  rep(0, J)
  Wcons0[as.numeric(names(table(pref_model)))] <- Wcons
  Wcons <- Wcons0}
(Wcons <- matrix(Wcons))


# plot weights data -------------------------------------------------------

# global & consensus weights dynamics
Wglobal_dyn <- apply(Wtemp, 3, colMeans)
Wcons_dyn <- matrix(0,tmax,J)
Wcons_dyn[1,] <- rep(1/J,J)
for (t in 1:(tmax-1)) {
  pref_model_t <- apply(Wtemp[,,t], 1, which.max)
  Wcons_dyn[t+1,sort(unique(pref_model_t))] <- (table(pref_model_t)/sum(table(pref_model_t)))
}
mod_labels <- apply(par_grid, 1, function(row) paste0("(", row[1], ", ", row[2], ")"))
rownames(Wglobal_dyn) <- mod_labels
colnames(Wcons_dyn)   <- mod_labels
df_global <- as.data.frame(Wglobal_dyn) %>%
  mutate(Component = mod_labels) %>%
  pivot_longer(-Component, names_to = "Time", values_to = "Weight") %>%
  mutate(Type = "W global")
df_cons <- as.data.frame(t(Wcons_dyn)) %>%
  mutate(Component = mod_labels) %>%
  pivot_longer(-Component, names_to = "Time", values_to = "Weight") %>%
  mutate(Type = "W consensus")
df_all <- bind_rows(df_global, df_cons) %>%
  mutate(Time = as.numeric(gsub("V", "", Time)))  # assuming default colnames are V1, V2, ...
pal <- brewer.pal(J, "Paired")

# Plot
weigth_dyn_plot <- ggplot(df_all, aes(x = Time, y = Weight, color = Component)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Type, ncol = 2, scales = "free_y") +
  scale_color_manual(values = pal) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Time", y = "Weight", color = "Candidate") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    # legend.position = "right",
    # legend.title = element_text(face = "bold", size = 14),
    # legend.text = element_text(face = "bold", size = 12),
    # legend.key = element_rect(fill = "white", color = NA),
    # legend.background = element_rect(fill = "gray95", color = "gray80"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )
print(weigth_dyn_plot)


# average estimated parameters over time
par_dyn <- sapply(2:tmax, function(t) matrix(colMeans(Wtemp[,,t]), nrow = 1) %*% par_grid)
par_df <- data.frame(
  Time = 2:tmax,
  alpha = par_dyn[1, ],
  phi = par_dyn[2, ]
) %>%
  pivot_longer(cols = c("alpha", "phi"), names_to = "Parameter", values_to = "Value") %>%
  mutate(
    TrueValue = case_when(
      Parameter == "alpha" ~ tau,
      Parameter == "phi" ~ phi
    )
  )

label_map <- c(alpha = bquote(alpha), phi = bquote(phi))
ylim_list <- list("alpha" = c(0.725, 0.825), "phi" = c(3.75, 5))
limit_df <- do.call(rbind, lapply(names(ylim_list), function(param) {
  data.frame(
    Parameter = param,
    Time = range(par_df$Time),
    Value = ylim_list[[param]]
  )
}))

# Plot
par_dyn_plot <- ggplot(par_df, aes(x = Time, y = Value)) +
  geom_line(aes(color = Parameter), linewidth = 1.2) +
  geom_hline(aes(yintercept = TrueValue), color = "#d7191c", linetype = "dashed", linewidth = 1) +
  geom_blank(data = limit_df) +  # 👈 Force facet limits
  facet_wrap(
    ~Parameter,
    scales = "free_y",
    labeller = as_labeller(label_map, default = label_parsed)
  ) +
  scale_color_manual(values = c("alpha" = "#6f59b4", "phi" = "#fdae61")) +
  labs(x = "Time", y = NULL) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(par_dyn_plot)

# global & consensus weights spatial distribution (voronoi tassellation)
pref_model <- apply(Wi, 1, which.max)

# compute Voronoi tessellation (deldir)
tess <- deldir(crd[,1], crd[,2], rw = c(-0.01, 1.01, -0.01, 1.01))
tiles <- tile.list(tess)
vor_df <- do.call(rbind, lapply(seq_along(tiles), function(i) {
  data.frame(
    x = tiles[[i]]$x,
    y = tiles[[i]]$y,
    id = i,
    pref_model = pref_model[i]
  )
}))

# voronoi plot with ggplot
vor_df$Model <- factor(vor_df$pref_model, levels = seq_along(mod_labels), labels = mod_labels)
points_df <- data.frame(x = crd[,1], y = crd[,2], Model = factor(pref_model, levels = seq_along(mod_labels), labels = mod_labels))
voronoi_plot <- ggplot() +
  geom_polygon(data = vor_df, aes(x = x, y = y, group = id, fill = Model), color = "black", alpha = 0.6) +
  geom_point(data = points_df, aes(x = x, y = y, color = Model), size = 2) +
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
  scale_fill_manual(values = pal, name = "Candidate") +
  scale_color_manual(values = pal, guide = "none") +
  labs(title = "", x = "Easting", y = "Northing") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    legend.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid = element_blank()
  )

# barplot
bar_data <- data.frame(
  Model = factor(mod_labels, levels = mod_labels),
  MeanWeight = colMeans(Wi)
)

bar_plot <- ggplot(bar_data, aes(x = Model, y = MeanWeight, fill = Model)) +
  geom_col() +
  scale_fill_manual(values = pal, name = "Candidate") +
  labs(title = "", x = "Model", y = "Average Weight") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    # legend.title = element_text(face = "bold", size = 14),
    # legend.text = element_text(face = "bold", size = 12),
    # legend.key = element_rect(fill = "white", color = NA),
    # legend.background = element_rect(fill = "gray95", color = "gray80"),
    legend.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

# combine
# combined_plot <- voronoi_plot + bar_plot + plot_layout(widths = c(2, 1)) & theme(legend.position = "right")

combined_plot <- (bar_plot) + (voronoi_plot + theme(legend.position = "right")) +
  plot_layout(widths = c(1, 2))

print(combined_plot)

# save plots --------------------------------------------------------------

save_plots_to_png <- function(
    plot_list,
    prefix = "plot",
    folder = ".",
    height = 3000,
    res = 300,
    ratio = 16/9,
    pointsize = 12
) {
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  
  for (i in seq_along(plot_list)) {
    # Sanitize plot name
    name_part <- if (!is.null(names(plot_list)) && names(plot_list)[i] != "") {
      gsub("[^A-Za-z0-9_-]", "_", names(plot_list)[i])
    } else {
      paste0("plot_", i)
    }
    
    fname <- file.path(folder, paste0(prefix, "_", name_part, ".png"))
    
    width = height*(ratio)
    png(filename = fname, width = width, height = height, res = res, pointsize = pointsize)
    print(plot_list[[i]])
    dev.off()
    
    message("✅ Saved PNG: ", fname)
  }
}


save_plots_to_png(list("weights_dynamic" = weigth_dyn_plot),
                  prefix = "plot", folder = "plots", height = 1250, ratio = 24/7)
save_plots_to_png(list("par_dynamic" = par_dyn_plot),
                  prefix = "plot", folder = "plots", height = 1250, ratio = 24/8)
save_plots_to_png(list( "weight_distr" = combined_plot),
                  prefix = "plot", folder = "plots", height = 1500, ratio = 24/12)


