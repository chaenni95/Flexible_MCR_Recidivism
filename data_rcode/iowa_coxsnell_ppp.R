library(rstan)
library(survival)
library(ggplot2)

iowa_rec18_dum <- read.csv("data_rcode/iowa18.csv")

llgfg.fit <- readRDS("results/data_result/llg4.rds")  # final model LLGFG

# x is for the survival time variable
# z is for the cure rate variable
iowa18_stan <- list(N = nrow(iowa_rec18_dum), t = iowa_rec18_dum$time_to_recurrence,
                    is_censored = iowa_rec18_dum$is_censored,
                    cured = iowa_rec18_dum$cured, x = cbind(1, iowa_rec18_dum[, -c(3,4,5,6,10, 13, 21, 22, 31, 32:38, 41:43)]),
                    z = iowa_rec18_dum[, -c(3,4,5,6,10, 13, 21, 22, 31, 32:43)], K = 22, M = 25)


X <- as.matrix(cbind(1, iowa_rec18_dum[, -c(3,4,5,6,10, 13, 21, 22, 31, 32:38, 41:43)]))
Z <- as.matrix(iowa_rec18_dum[, -c(3,4,5,6,10, 13, 21, 22, 31, 32:43)])

time <- iowa_rec18_dum$time_to_recurrence
status <- iowa_rec18_dum$is_censored


post.fit.llgfg <- rstan::extract(llgfg.fit)


### Cox-snell residual plot with final model #########################################################################
plot_cs_envelope_single <- function(fit,
                                    time,
                                    status,
                                    X,
                                    Z,
                                    dist = "weibull",
                                    link = "logit",
                                    tag = NULL,
                                    save_path = NULL,
                                    envelope = TRUE) {
  if (is.null(tag)) tag <- paste(dist, link, sep = "_")
  
  fit_post <- rstan::extract(fit)
  n <- length(time)
  L <- length(fit_post$sigma)
  S_mix_mat <- matrix(NA, n, L)
  
  for (l in 1:L) {
    betaU_l <- fit_post$betaU[l, ]
    betaC_l <- fit_post$betaC[l, ]
    sigma_l <- fit_post$sigma[l]
    
    delta_l <- if (!is.null(fit_post$delta)) fit_post$delta[l] else 1
    lambda_link_l <- if (!is.null(fit_post$lambda)) fit_post$lambda[l] else 1
    
    muU_l <- as.numeric(X %*% betaU_l)
    zbeta <- Z %*% betaC_l
    
    p_cure <- switch(link,
                     "logit" = 1 / (1 + exp(-zbeta)),
                     "slogit" = (1 / (1 + exp(-zbeta)))^delta_l,
                     "rplogit" = 1 - (1 / (1 + exp(zbeta)))^delta_l,
                     "fglogit" = (1 - (1 / (1 + exp(zbeta)))^lambda_link_l)^delta_l,
                     stop("Unsupported link"))
    
    Su_l <- switch(dist,
                   "weibull" = {
                     alpha_l <- 1 / sigma_l
                     lambda_l <- exp(muU_l)
                     exp(-(time / lambda_l)^alpha_l)
                   },
                   "lognormal" = {
                     lambda_l <- muU_l
                     1 - plnorm(time, meanlog = lambda_l, sdlog = sigma_l)
                   },
                   "loglogistic" = {
                     shape <- 1 / sigma_l
                     scale <- exp(muU_l)
                     1 / (1 + (time / scale)^shape)
                   },
                   stop("Unsupported distribution"))
    
    S_mix_mat[, l] <- p_cure + (1 - p_cure) * Su_l
  }
  
  r_cs_pp <- -log(rowMeans(S_mix_mat))
  km_fit <- survfit(Surv(r_cs_pp, status) ~ 1)
  km_df <- data.frame(time = km_fit$time, surv = km_fit$surv)
  
  grid_times <- seq(0, max(r_cs_pp), length.out = 100)
  exp_curve <- data.frame(time = grid_times, surv = exp(-grid_times))
  
  library(ggplot2)
  
  # Core KM data
  km_df <- data.frame(
    time = km_fit$time,
    surv = km_fit$surv,
    lower = km_fit$lower,
    upper = km_fit$upper
  )
  grid_times <- seq(0, max(r_cs_pp), length.out = 100)
  exp_curve <- data.frame(time = grid_times, surv = exp(-grid_times))
  
  # Envelope (if requested)
  if (envelope) {
    km_surv_matrix <- matrix(NA, nrow = length(grid_times), ncol = L)
    for (l in 1:L) {
      r_cs_l <- -log(S_mix_mat[, l])
      km_l <- survfit(Surv(r_cs_l, status) ~ 1)
      interp <- approx(km_l$time, km_l$surv, xout = grid_times, method = "linear", rule = 2)$y
      km_surv_matrix[, l] <- interp
    }
    
    lower <- apply(km_surv_matrix, 1, quantile, probs = 0.025, na.rm = TRUE)
    upper <- apply(km_surv_matrix, 1, quantile, probs = 0.975, na.rm = TRUE)
    envelope_df <- data.frame(time = grid_times, lower = lower, upper = upper)
  } else {
    envelope_df <- NULL
  }
  
  # Construct base ggplot
  p <- ggplot() +
    # Posterior envelope
    {if (!is.null(envelope_df)) geom_ribbon(data = envelope_df,
                                            aes(x = time, ymin = lower, ymax = upper, fill = "95% Envelope"),
                                            alpha = 0.5)} +
    # KM confidence interval (dotted lines)
    geom_step(data = km_df, aes(x = time, y = lower), color = "black", linetype = "dotted") +
    geom_step(data = km_df, aes(x = time, y = upper), color = "black", linetype = "dotted") +
    # KM curve
    geom_step(data = km_df, aes(x = time, y = surv, color = "KM of Residuals"), linewidth = 1) +
    # Exp(1) curve
    geom_line(data = exp_curve, aes(x = time, y = surv, color = "Exponential(1)"),
              linetype = "dashed", linewidth = 1) +
    labs(x = "Cox–Snell Residual", y = "Survival") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(name = NULL,
                       values = c("KM of Residuals" = "black",
                                  "Exponential(1)" = "red")) +
    scale_fill_manual(name = NULL,
                      values = c("95% Envelope" = "gray80")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top")
  
  
  # Save if needed
  if (!is.null(save_path)) {
    ggsave(filename = file.path(save_path, paste0("coxsnell_", tag, ".pdf")),
           plot = p, width = 7, height = 5)
  }
  return(list(S_mix_mat = S_mix_mat, p = p))
}




## Draw cox-snell plot for the final model
plot.llg4 <- plot_cs_envelope_single(
  fit = llgfg.fit,
  time = iowa_rec18_dum$time_to_recurrence,
  status = iowa_rec18_dum$is_censored,
  X = X,
  Z = Z,
  dist = "loglogistic",
  link = "fglogit",
  envelope = TRUE, 
  save_path = "results/data_result/", 
  tag = "loglogistic_fglogit"
)




### Posterior predictive p-value function #########################################################################

compute_pp_value <- function(S_mix_mat, status, grid_times = seq(0, 2, length.out = 100)) {
  L <- ncol(S_mix_mat)
  T_rep <- numeric(L)
  
  # Compute T_obs
  r_obs <- -log(rowMeans(S_mix_mat))
  km_obs <- survfit(Surv(r_obs, status) ~ 1)
  s_obs <- approx(km_obs$time, km_obs$surv, xout = grid_times, method = "constant", rule = 2)$y
  T_obs <- sum((s_obs - exp(-grid_times))^2)
  
  # Replicates
  for (l in 1:L) {
    r_l <- -log(S_mix_mat[, l])
    km_l <- survfit(Surv(r_l, status) ~ 1)
    s_l <- approx(km_l$time, km_l$surv, xout = grid_times, method = "constant", rule = 2)$y
    T_rep[l] <- sum((s_l - exp(-grid_times))^2)
  }
  
  ppp <- mean(T_rep >= T_obs)
  return(list(ppp = ppp, T_obs = T_obs, T_rep = T_rep))
}


result.llg4 <- compute_pp_value(plot.llg4$S_mix_mat, iowa_rec18_dum$is_censored)
result.llg4$ppp  # Posterior predictive p-value
0.5855