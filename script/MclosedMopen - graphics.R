# Packages ----------------------------------------------------------------

library(ggplot2)
library(ggridges)
library(dplyr)
library(patchwork)
library(tidyr)
library(xtable)


# Sigma -------------------------------------------------------------------


method_names <- c("True", "DynBPS-C", "DynBPS-O")
method_colors <- c(
  "True"     = "#fdd835",   # Strong red
  "DynBPS-C" = "#1976d2",   # Strong blue
  "DynBPS-O" = "#2ecc71"    # Emerald
)
# "#d7191c"

true_cov <- Sigma

results_all <- data.frame()
for (m in seq_along(method_names)) {
  
  method <- method_names[m]
  # posterior_samples <- posterior_list[[m]]  # Assuming a list of arrays per method
  
  sigma_post <- sapply(out_list, function(b){ b[[m]]$`Sigma posterior` }, simplify = "array")
  posterior_samples <- sigma_post
  
  method_results <- data.frame()
  
  for (i in 1:q) {
    for (j in i:q) {
      est_means <- numeric(BB)
      est_sds <- numeric(BB)
      coverages <- logical(BB)
      
      for (r in 1:BB) {
        samples <- posterior_samples[i, j, , r]
        est_means[r] <- mean(samples)
        est_sds[r] <- sd(samples)
        ci <- quantile(samples, probs = c(0.025, 0.975))
        coverages[r] <- (true_cov[i, j] >= ci[1] && true_cov[i, j] <= ci[2])
      }
      
      method_results <- rbind(method_results, data.frame(
        Method = method,
        i = i,
        j = j,
        Entry = paste0("(", i, ",", j, ")"),
        Bias = mean(est_means - true_cov[i, j]),
        Coverage = mean(coverages),
        SD = mean(est_sds)
      ))
    }
  }
  
  results_all <- rbind(results_all, method_results)
}

entries <- c("(1,1)", "(1,2)", "(2,2)", "(1,3)", "(2,3)", "(3,3)")
greek_labels <- c("Sigma[1,1]", "Sigma[1,2]", "Sigma[2,2]", "Sigma[1,3]", "Sigma[2,3]", "Sigma[3,3]")
param_labels <- c(
  expression(Sigma[1*","*1]), expression(Sigma[1*","*2]), expression(Sigma[2*","*2]),
  expression(Sigma[1*","*3]), expression(Sigma[2*","*3]), expression(Sigma[3*","*3])
)
names(param_labels) <- greek_labels
results_all$Label <- factor(results_all$Entry, levels = entries, labels = greek_labels)

sigma_bias <- ggplot(results_all, aes(x = Label, y = Bias, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.35) +
  coord_cartesian(ylim = c(-.35, .35)) +
  scale_fill_manual(values = method_colors) +
  scale_x_discrete(labels = param_labels) +
  labs(x = "", y = "Avg. Empirical Bias") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    # legend.position = "right"
    legend.position = "none"
  ) 

sigma_cover <- ggplot(results_all, aes(x = Label, y = Coverage, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.35) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", linewidth = 1.2) +
  coord_cartesian(ylim = c(0, 1.2)) +
  scale_fill_manual(values = method_colors) +
  scale_x_discrete(labels = param_labels) +
  labs(x = "", y = "Avg. Coverage") +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      plot.title = element_text(size = 16, face = "bold"),
      # legend.position = "right"
      legend.position = "none"
    )

sigma_std <- ggplot(results_all, aes(x = Label, y = SD, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.35) +
  coord_cartesian(ylim = c(0, .3)) +
  scale_fill_manual(values = method_colors) +
  scale_x_discrete(labels = param_labels) +
  labs(x = "", y = "Avg. Std. Dev") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    # legend.position = "right"
    legend.position = "none"
  )


frobenius_df <- data.frame()

for (m in seq_along(method_names)) {
  
  method <- method_names[m]
  # posterior_samples <- posterior_list[[m]]  # Assuming a list of arrays per method
  
  sigma_post <- sapply(out_list, function(b){ b[[m]]$`Sigma posterior` }, simplify = "array")
  posterior_samples <- sigma_post
  
  norms <- numeric(BB)
  for (r in 1:BB) {
    post_mean <- apply(posterior_samples[, , , r], c(1,2), mean)
    diff <- post_mean - true_cov
    norms[r] <- sqrt(sum(diff^2))
  }
  
  frobenius_df <- rbind(frobenius_df, data.frame(
    Replication = 1:BB,
    FrobeniusNorm = norms,
    Method = method
  ))
}

sigma_frob <- ggplot(frobenius_df, aes(x = Method, y = FrobeniusNorm, fill = Method)) +
  geom_violin(trim = F) +
  geom_boxplot(width = 0.1, outliers = F) +
  coord_cartesian(ylim = c(0, .75)) +
  scale_fill_manual(values = method_colors) +
  labs(y = bquote("Frobenius norm : || " * bar(Sigma)  - Sigma[true] * " || "[F]), x = "") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )


sigma_plot <- ((sigma_bias / sigma_cover / sigma_std + plot_layout(heights = c(1, 1, 1)))
               | sigma_frob) + plot_layout(widths = c(1.5, 1))
print(sigma_plot)


# Theta -------------------------------------------------------------------

method_names <- c("True", "DynBPS-C", "DynBPS-O")
method_colors <- c(
  "True"     = "#fdd835",   # Strong red
  "DynBPS-C" = "#1976d2",   # Strong blue
  "DynBPS-O" = "#2ecc71"    # Emerald
)

theta_true <- theta[1:p,,]
results_theta <- data.frame()

for (m in seq_along(method_names)) {
  
  method <- method_names[m]
  posterior <- sapply(out_list, function(b){ b[[m]]$`theta posterior` }, simplify = "array")
  
  for (i in 1:p) {
    for (j in 1:q) {
      for (t in 1:tmax) {
        est_means <- numeric(BB)
        est_sds <- numeric(BB)
        coverages <- logical(BB)
        
        for (r in 1:BB) {
          samples <- posterior[i, j, , t, r]
          true_val <- theta_true[i, j, t]
          est_means[r] <- mean(samples)
          est_sds[r] <- sd(samples)
          ci <- quantile(samples, probs = c(0.025, 0.975))
          coverages[r] <- (true_val >= ci[1] && true_val <= ci[2])
        }
        
        results_theta <- rbind(results_theta, data.frame(
          Method = method,
          i = i,
          j = j,
          Entry = paste0("(", i, ",", j, ")"),
          Time = t,
          Bias = mean(est_means - true_val),
          SD = mean(est_sds),
          Coverage = mean(coverages)
        ))
      }
    }
  }
}

entries <- c("(1,1)", "(1,2)", "(1,3)", "(2,1)", "(2,2)", "(2,3)")
param_labels <- c(
  expression(Theta[1*","*1]), expression(Theta[1*","*2]), expression(Theta[1*","*3]),
  expression(Theta[2*","*1]), expression(Theta[2*","*2]), expression(Theta[2*","*3])
)
results_theta$Label <- factor(results_theta$Entry, levels = entries, labels = param_labels)


theta_bias <- ggplot(results_theta, aes(x = factor(Time), y = Bias, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ Label, labeller = label_parsed) +
  coord_cartesian(ylim = c(-.5, .5)) +
  scale_fill_manual(values = method_colors) +
  labs(x = "Time", y = "Avg. Empirical Bias") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    # legend.position = "right"
    legend.position = "none"
  )

theta_cover <- ggplot(results_theta, aes(x = factor(Time), y = Coverage, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", linewidth = 1.2) +
  facet_wrap(~ Label, labeller = label_parsed) +
  coord_cartesian(ylim = c(0, 1.2)) +
  scale_fill_manual(values = method_colors) +
  labs(x = "Time", y = "Avg. Coverage") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    # legend.position = "right"
    legend.position = "none"
  )


theta_std <- ggplot(results_theta, aes(x = factor(Time), y = SD, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ Label, labeller = label_parsed) +
  coord_cartesian(ylim = c(0, .5)) +
  scale_fill_manual(values = method_colors) +
  labs(x = "Time", y = "Avg. Std. Dev") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    # legend.position = "right"
    legend.position = "none"
  )


frob_theta <- data.frame()

for (m in seq_along(method_names)) {
  
  method <- method_names[m]
  posterior <- sapply(out_list, function(b){ b[[m]]$`theta posterior` }, simplify = "array")
  
  for (r in 1:BB) {
    for (t in 1:tmax) {
      post_mean <- apply(posterior[1:p, , , t, r], c(1,2), mean)
      true_mat <- theta_true[, , t]
      norm <- sqrt(sum((post_mean - true_mat)^2))
      
      frob_theta <- rbind(frob_theta, data.frame(
        Method = method,
        Replication = r,
        Time = t,
        FrobeniusNorm = norm
      ))
    }
  }
}

medians_df <- frob_theta %>%
  group_by(Method, Time) %>%
  summarize(median_frob = median(FrobeniusNorm), .groups = "drop") %>%
  mutate(Time_factor = factor(Time))

theta_frob <- ggplot(frob_theta, aes(x = factor(Time), y = FrobeniusNorm, fill = Method)) +
  geom_violin(trim = T, width = 1) +
  geom_boxplot(width = 0.2, outliers = F, position = position_dodge(0.9)) +
  geom_point(
    data = medians_df,
    aes(x = Time_factor, y = median_frob, color = Method),
    size = 1.5,
    position = position_dodge(width = 0.5)
  ) +
  geom_line(
    data = medians_df,
    aes(x = Time_factor, y = median_frob, color = Method, group = Method),
    linewidth = 0.8,
    position = position_dodge(width = 0.5)
  ) +
  coord_cartesian(ylim = c(0, 1.5)) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  labs(x = "Time", y = expression("Frobenius norm : ||" * Theta - Theta[true] * "||"[F])) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )
  

theta_plot <- (theta_bias / theta_cover) | (theta_std / theta_frob)
print(theta_plot)

# Omega -------------------------------------------------------------------

method_names <- c("True", "DynBPS-C", "DynBPS-O")
method_colors <- c(
  "True"     = "#fdd835",   # Strong red
  "DynBPS-C" = "#1976d2",   # Strong blue
  "DynBPS-O" = "#2ecc71"    # Emerald
)

omega_true <- theta
results_omega <- data.frame()

for (m in seq_along(method_names)) {
  
  method <- method_names[m]
  posterior <- sapply(out_list, function(b){ b[[m]]$`theta posterior` }, simplify = "array")
  
  for (i in (p+1):(p+n)) {
    for (j in 1:q) {
      for (t in 1:tmax) {
        est_means <- numeric(BB)
        est_sds <- numeric(BB)
        coverages <- logical(BB)
        
        for (r in 1:BB) {
          samples <- posterior[i, j, , t, r]
          true_val <- omega_true[i, j, t]
          est_means[r] <- mean(samples)
          est_sds[r] <- sd(samples)
          ci <- quantile(samples, probs = c(0.025, 0.975))
          coverages[r] <- (true_val >= ci[1] && true_val <= ci[2])
        }
        
        results_omega <- rbind(results_omega, data.frame(
          Method = method,
          i = i,
          j = j,
          Entry = paste0("(", i, ",", j, ")"),
          Time = t,
          Bias = mean(est_means - true_val),
          SD = mean(est_sds),
          Coverage = mean(coverages)
        ))
      }
    }
  }
}


omega_bias <- ggplot(results_omega, aes(x = Method, y = Bias, fill = Method)) +
  geom_violin(trim = T) +
  geom_boxplot(width = 0.2, outliers = F, position = position_dodge(0.9)) +
  facet_wrap(~ j, nrow = 1, labeller = label_bquote(Omega[.(j)])) +
  scale_fill_manual(values = method_colors) +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "Time", y = "Avg. Empirical Bias") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    # legend.position = "right"
    legend.position = "none"
  )

omega_cover <- ggplot(results_omega, aes(x = Method, y = Coverage, fill = Method)) +
  geom_violin(trim = T, bounds = c(0.85, 1)) +
  geom_boxplot(width = 0.2, outliers = F, position = position_dodge(0.9)) +
  facet_wrap(~ j, nrow = 1, labeller = label_bquote(Omega[.(j)])) +
  scale_fill_manual(values = method_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Time", y = "Avg. Coverage") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    # legend.position = "right"
    legend.position = "none"
  )

omega_std <- ggplot(results_omega, aes(x = Method, y = SD, fill = Method)) +
  geom_violin(trim = T) +
  geom_boxplot(width = 0.2, outliers = F, position = position_dodge(0.9)) +
  facet_wrap(~ j, nrow = 1, labeller = label_bquote(Omega[.(j)])) +
  scale_fill_manual(values = method_colors) +
  labs(x = "Time", y = "Avg. Std. Dev.") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    # legend.position = "right"
    legend.position = "none"
  )


frob_omega <- data.frame()

for (m in seq_along(method_names)) {
  
  method <- method_names[m]
  posterior <- sapply(out_list, function(b){ b[[m]]$`theta posterior` }, simplify = "array")
  
  for (r in 1:BB) {
    for (t in 1:tmax) {
      post_mean <- apply(posterior[-(1:p), , , t, r], c(1,2), mean)
      true_mat <- omega_true[(p+1):(p+n), , t]
      norm <- sqrt(sum((post_mean - true_mat)^2))
      
      frob_omega <- rbind(frob_omega, data.frame(
        Method = method,
        Replication = r,
        Time = t,
        FrobeniusNorm = norm
      ))
    }
  }
}

# Compute mean Frobenius norms by Method and Time
medians_df <- frob_omega %>%
  group_by(Method, Time) %>%
  summarize(median_frob = median(FrobeniusNorm), .groups = "drop") %>%
  mutate(Time_factor = factor(Time))

omega_frob <- ggplot(frob_omega, aes(x = factor(Time), y = FrobeniusNorm, fill = Method)) +
  geom_violin(trim = T, width = 1) +
  geom_boxplot(width = 0.2, outliers = F, position = position_dodge(0.9)) +
  geom_point(
    data = medians_df,
    aes(x = Time_factor, y = median_frob, color = Method),
    size = 1.5,
    position = position_dodge(width = 0.5)
  ) +
  geom_line(
    data = medians_df,
    aes(x = Time_factor, y = median_frob, color = Method, group = Method),
    linewidth = 0.8,
    position = position_dodge(width = 0.5)
  ) +
  coord_cartesian(ylim = c(6, 12)) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  labs(x = "Time", y = expression("Frobenius norm : ||" * Omega - Omega[true] * "||"[F])) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )


omega_plot <- (omega_bias / omega_cover) | (omega_std / omega_frob)
print(omega_plot)

# library(knitr)
# library(kableExtra)
# agg_table %>%
#   kable("latex", digits = 3, booktabs = TRUE,
#         caption = "Summary of Predictive Metrics for $\\Omega$ by Method and Variable $j$") %>%
#   kable_styling(latex_options = c("hold_position", "striped"))



# Predictions -------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(knitr)

method_names <- c("True", "DynBPS-C", "DynBPS-O")
method_colors <- c(
  "True"     = "#fdd835",   # Strong red
  "DynBPS-C" = "#1976d2",   # Strong blue
  "DynBPS-O" = "#2ecc71"    # Emerald
)

test <- sapply(out_list, function(b) b[[4]], simplify = "array")

evaluate_all_methods <- function(true_data, forecast_only = FALSE, out_list, TMAX = NULL) {
  metrics_df <- list()
  
  for (m in seq_along(method_names)) {
    method <- method_names[m]
    pred_array <- sapply(out_list, function(b) b[[m]]$prediction, simplify = "array")
    dims <- dim(pred_array)
    n <- dims[1]; q <- dims[2]; S <- dims[3]; TT <- dims[4]; R <- dims[5]
    
    time_range <- if (forecast_only) (TMAX + 1):TT else 1:TT
    
    for (r in 1:R) {
      for (t in time_range) {
        for (j in 1:q) {
          mspe_vec <- numeric(n)
          bias_vec <- numeric(n)
          pi_width_vec <- numeric(n)
          var_vec <- numeric(n)
          
          for (i in 1:n) {
            pred_samples <- pred_array[i, j, , t, r]
            true_val <- true_data[i, j, t, r]
            
            mean_pred <- mean(pred_samples)
            ci <- quantile(pred_samples, c(0.025, 0.975))
            
            mspe_vec[i] <- (mean_pred - true_val)^2
            bias_vec[i] <- abs(mean_pred - true_val)
            pi_width_vec[i] <- ci[2] - ci[1]
            var_vec[i] <- var(pred_samples)
          }
          
          metrics_df[[length(metrics_df) + 1]] <- data.frame(
            Method = method,
            Time = t,
            Replication = r,
            Variable = as.factor(j),
            MSPE = mean(mspe_vec),
            Bias = mean(bias_vec),
            PI_Width = mean(pi_width_vec),
            Pred_Var = mean(var_vec)
          )
        }
      }
    }
  }
  
  bind_rows(metrics_df) %>%
    pivot_longer(cols = c("MSPE", "Bias", "PI_Width", "Pred_Var"),
                 names_to = "Metric", values_to = "Value")
}

plot_box <- function(df, tmax) {
  df <- df %>%
    mutate(
      Metric = case_when(
        Metric %in% c("MSPE", "mspe") ~ "MSPE",
        Metric %in% c("Bias", "bias") ~ "Abs. Bias",
        Metric %in% c("PI_Width", "Pred. Int. Width", "PI.Width") ~ "Pred. Int. Width",
        Metric %in% c("Pred_Var", "Pred. Var.") ~ "Pred. Variance",
        TRUE ~ as.character(Metric)
      ),
      Variable = factor(Variable, labels = paste0("Y[", levels(factor(Variable)), "]")),
      Horizon = as.factor(as.integer(Time) - tmax)
    )
  
  ggplot(df, aes(x = Horizon, y = Value, fill = Method)) +
    geom_boxplot(width = 0.3, outliers = F) +
    facet_grid(Metric ~ Variable, scales = "free_y", 
               labeller = labeller(Variable = label_parsed)) +
    scale_fill_manual(values = method_colors) +
    labs(x = "Forecast Horizon", y = NULL) + 
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "right"
    )
}

summary_table <- function(df) {
  df %>%
    group_by(Method, Metric) %>%
    summarise(Mean = mean(Value), SD = sd(Value), .groups = "drop") %>%
    mutate(across(where(is.numeric), ~round(.x, 4))) %>%
    pivot_wider(names_from = Metric, values_from = c(Mean, SD)) %>%
    kable(format = "latex", booktabs = TRUE, caption = "Summary of Predictive Metrics")
}

# Usage
metrics_long <- evaluate_all_methods(true_data = test[1:n,,,], forecast_only = T, out_list = out_list, TMAX = tmax)

# Visualize
pred_plot <- plot_box(metrics_long, tmax = tmax)
print(pred_plot)

# Summary
summary_table(metrics_long)


# Save plots --------------------------------------------------------------

save_plots_to_png <- function(
    plot_list,
    prefix = "plot",
    folder = ".",
    height = 3000,
    res = 300,
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
    
    width = height*(16/9)
    png(filename = fname, width = width, height = height, res = res, pointsize = pointsize)
    print(plot_list[[i]])
    dev.off()
    
    message("✅ Saved PNG: ", fname)
  }
}


save_plots_to_png(list("sigma" = sigma_plot, "theta" = theta_plot, "omega" = omega_plot, "pred" = pred_plot),
                  prefix = "plot", folder = "plots")


