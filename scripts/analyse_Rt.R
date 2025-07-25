analyse_Rt <- function(incidence, start_date, end_date, window_width, mean_si, std_si) {
  
  # Set the ggplot theme
  plot_theme <- theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 15)
    )
  
  incidence_subset <- incidence |> filter(Date >= start_date & Date <= end_date)
  n_obs <- nrow(incidence_subset) # How many observations

  # Calculate the starting and end points of the time windows
  t_starts <- 2:(n_obs - window_width + 1)
  t_ends <- t_starts + window_width - 1
  # Discretization taken from EpiEstim
  si_distr <- discr_si(seq(0, nrow(incidence) - 1), mean_si, std_si)
  # Compute lambda
  lambda <- vector(mode = "numeric", length = nrow(incidence))
  lambda[1] <- NA
  for (t in seq(2, nrow(incidence))) {
    lambda[t] <- sum(
      si_distr[seq_len(t)] * incidence$Cases[seq(t, 1)]
    )
  }
  lambda_subset <- lambda[incidence$Date >= start_date & incidence$Date <= end_date]

  pairs <- cbind(lambda_subset, incidence_subset$Cases)
  colnames(pairs) <- c("lambda", "Cases")

  models_pois <- models_qpois <- models_nbin_Q <- models_nbin_L <-
    vector(mode = "list", length = length(t_starts))
  model_mats <- sapply(
    t_starts,
    function(x) {
      pairs[0:(window_width - 1) + x, ]
    },
    simplify = FALSE
  ) %>%
    lapply(as.data.frame)

  for (k in 1:length(t_starts)) {
    # Poisson model
    models_pois[[k]] <- glm(
      data = model_mats[[k]],
      Cases ~ lambda - 1,
      family = poisson(link = "identity")
    )

    # Quasipoisson model
    models_qpois[[k]] <- glm(
      data = model_mats[[k]],
      Cases ~ lambda - 1,
      family = quasipoisson(link = "identity")
    )

    # NegBin-Q model
    models_nbin_Q[[k]] <- gamlss(
      data = model_mats[[k]],
      formula = Cases ~ lambda - 1,
      family = NBI(mu.link = "identity", sigma.link = "log"),
      control = gamlss.control(trace = FALSE)
    )

    # NegBin-L model
    models_nbin_L[[k]] <- gamlss(
      data = model_mats[[k]],
      formula = Cases ~ lambda - 1,
      family = NBII(mu.link = "identity", sigma.link = "log"),
      control = gamlss.control(trace = FALSE)
    )
  }

  # Extract the coefficient
  R_hat <- list()
  R_hat$pois <- models_pois %>%
    lapply(coefficients) %>%
    unlist()
  R_hat$qpois <- models_qpois %>%
    lapply(coefficients) %>%
    unlist()
  R_hat$nbin_Q <- models_nbin_Q %>%
    lapply(coefficients) %>%
    unlist()
  R_hat$nbin_L <- models_nbin_L %>%
    lapply(coefficients) %>%
    unlist()
  # NegBin-Q model, estimates based on the approximate formulas
  R_hat$nbin_Q_approx <- model_mats |> 
    lapply(function (x) {sum(x$Cases / x$lambda) / window_width}) |> unlist()

  # Extract the dispersion parameter estimates
  disp <- list()
  disp$qpois <- models_qpois %>%
    lapply(summary) %>%
    lapply(function(x) x$dispersion) %>%
    unlist()
  disp$nbin_Q <- models_nbin_Q %>%
    lapply(function(x) {
      exp(x$sigma.coefficients[1])
    }) %>%
    unlist() %>%
    unname()
  disp$nbin_L <- models_nbin_L %>%
    lapply(function(x) {
      exp(x$sigma.coefficients[1])
    }) %>%
    unlist() %>%
    unname()
  # NegBin-Q model, estimates based on the approximate formulas
  disp$nbin_Q_approx <- vapply(
    X = 1:length(t_starts),
    FUN.VALUE = 0,
    FUN = function (j) {
      num <- (model_mats[[j]]$Cases - R_hat$nbin_Q_approx[j] * model_mats[[j]]$lambda)^2
      denom <- (R_hat$nbin_Q_approx[j] * model_mats[[j]]$lambda)^2
      sum(num / denom) / (window_width - 1)
    }
  )
  
  # Extract the standard error
  R_hat_sd <- list()
  R_hat_sd$pois <- models_pois %>%
    lapply(summary) %>%
    lapply(function(x) x$coefficients[2]) %>%
    unlist()
  R_hat_sd$qpois <- models_qpois %>%
    lapply(summary) %>%
    lapply(function(x) x$coefficients[2]) %>%
    unlist()
  R_hat_sd$nbin_Q <- models_nbin_Q %>%
    lapply(summary) %>%
    lapply(function(x) x[3]) %>%
    unlist()
  R_hat_sd$nbin_L <- models_nbin_L %>%
    lapply(summary) %>%
    lapply(function(x) x[3]) %>%
    unlist()
  R_hat_sd$nbin_Q_approx <- sqrt(R_hat$nbin_Q_approx^2 * disp$nbin_Q_approx / window_width)

  # Get the AIC values
  AIC_vals <- list()
  AIC_vals$pois <- models_pois %>%
    lapply(function(x) x$aic) %>%
    unlist()
  AIC_vals$nbin_Q <- models_nbin_Q %>%
    lapply(function(x) x$aic) %>%
    unlist()
  AIC_vals$nbin_L <- models_nbin_L %>%
    lapply(function(x) x$aic) %>%
    unlist()

  # Create plots ----------------------------------------------------
  df_R_hat <- tibble(
    Date = rep(incidence_subset$Date[t_ends], 4),
    R = c(R_hat$pois, R_hat$qpois, R_hat$nbin_Q, R_hat$nbin_L),
    # Wald CI
    lwr = c(
      R_hat$pois - qnorm(0.975) * R_hat_sd$pois,
      R_hat$qpois - qnorm(0.975) * R_hat_sd$qpois,
      R_hat$nbin_Q - qnorm(0.975) * R_hat_sd$nbin_Q,
      R_hat$nbin_L - qnorm(0.975) * R_hat_sd$nbin_L
    ),
    upr = c(
      R_hat$pois + qnorm(0.975) * R_hat_sd$pois,
      R_hat$qpois + qnorm(0.975) * R_hat_sd$qpois,
      R_hat$nbin_Q + qnorm(0.975) * R_hat_sd$nbin_Q,
      R_hat$nbin_L + qnorm(0.975) * R_hat_sd$nbin_L
    ),
    Model = factor(rep(
      c("Poiss", "Q-Poiss", "NegBin-Q", "NegBin-L"),
      each = length(models_pois)
    ))
  )

  # Create plots ----------------------------------------------------
  # These colors are selected from the ggokabeito package designed for creating
  # color-blindness friendly charts
  model_colors <- c(
    "Poiss" = "#009E73", "Q-Poiss" = "#F0E442",
    "NegBin-L" = "#56B4E9", "NegBin-Q" = "#CC79A7"
  )
  
  # 1. Incidence
  p_incidence <- ggplot(
    incidence_subset,
    aes(x = Date, y = Cases)
  ) +
    geom_line(key_glyph = "timeseries", linewidth = 0.8) +
    scale_color_manual(values = "black", name = "") +
    labs(title = "Incidence") +
    plot_theme +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(0, max(incidence_subset$Cases)))

  # 2. Poisson vs. QP
  p_pois_vs_qpois <- ggplot(
    df_R_hat,
    aes(
      x = Date, y = R, ymin = lwr, ymax = upr, color = Model, fill = Model,
      alpha = Model
    )
  ) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_ribbon(color = NA) +
    scale_alpha_manual(
      values = c("Poiss" = 0.5, "Q-Poiss" = 0.5, "NegBin-L" = 0.0, "NegBin-Q" = 0),
      guide = guide_legend(override.aes = list(alpha = 0.5))
    ) +
    scale_color_manual(
      name = "Model",
      values = model_colors
    ) +
    scale_fill_manual(
      name = "Model",
      values = model_colors
    ) +
    labs(
      title = "Poisson vs. Quasi-Poisson",
      y = expression(hat(R))
    ) +
    coord_cartesian(
      ylim = c(0, 3), 
      xlim = range(incidence_subset$Date)
    ) +
    plot_theme

  # 3. NegBin-L vs. NegBin-Q
  p_nbin_L_vs_nbin_Q <- ggplot(
    df_R_hat,
    aes(
      x = Date, y = R, ymin = lwr, ymax = upr, color = Model, fill = Model,
      alpha = Model
    )
  ) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_ribbon(color = NA) +
    scale_alpha_manual(
      values = c("Poiss" = 0, "Q-Poiss" = 0, "NegBin-L" = 0.5, "NegBin-Q" = 0.5),
      guide = guide_legend(override.aes = list(alpha = 0.5))
    ) +
    scale_color_manual(
      name = "Model",
      values = model_colors
    ) +
    scale_fill_manual(
      name = "Model",
      values = model_colors
    ) +
    labs(
      title = "NegBin-L vs. NegBin-Q",
      y = expression(hat(R))
    ) +
    coord_cartesian(
      ylim = c(0, 3), 
      xlim = range(incidence_subset$Date)
    ) +
    plot_theme

  # 4. Overdisp Params. NegBin-Q plotted on using a second plot axis
  sec_axis_scale <- median(disp$nbin_Q / disp$nbin_L)
  disp$nbin_Q_transformed <- disp$nbin_Q / sec_axis_scale
  df_disp <- tibble(
    Date = rep(incidence_subset$Date[t_ends], 3),
    Dispersion = c(disp$qpois, disp$nbin_L, disp$nbin_Q_transformed),
    Model = factor(rep(c("Q-Poiss", "NegBin-L", "NegBin-Q"), each = length(disp$qpois)))
  )

  p_disp <- ggplot(
    df_disp,
    aes(x = Date, y = Dispersion, color = Model, group = Model)
  ) +
    geom_line(linewidth = 0.8) +
    scale_y_continuous(sec.axis = sec_axis(~ . * sec_axis_scale, name = "NegBin-Q")) +
    scale_color_manual(
      name = "Model",
      values = model_colors[-1]
    ) +
    labs(
      title = "Overdispersion Parameters",
      y = "Overdispersion"
    ) +
    coord_cartesian(
      ylim = c(0, max(df_disp$Dispersion)), 
      xlim = range(incidence_subset$Date)
    ) +
    plot_theme

  # 5. NegBin-L vs. QP (Extra plot, will probably go to the supplement)
  p_nbin_L_vs_qpois <- ggplot(
    df_R_hat,
    aes(
      x = Date, y = R, ymin = lwr, ymax = upr, color = Model, fill = Model,
      alpha = Model
    )
  ) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_ribbon(color = NA) +
    scale_alpha_manual(
      values = c("Poiss" = 0, "Q-Poiss" = 0.5, "NegBin-L" = 0.5, "NegBin-Q" = 0),
      breaks = c("Q-Poiss", "NegBin-L"),
      guide = guide_legend(override.aes = list(alpha = 0.5))
    ) +
    scale_color_manual(
      name = "Model",
      values = model_colors,
      breaks = c("Q-Poiss", "NegBin-L")
    ) +
    scale_fill_manual(
      name = "Model",
      values = model_colors,
      breaks = c("Q-Poiss", "NegBin-L")
    ) +
    labs(
      title = "NegBin-L vs. Quasi-Poisson",
      y = expression(hat(R))
    ) +
    coord_cartesian(
      ylim = c(0, 3), 
      xlim = range(incidence_subset$Date)
    ) +
    plot_theme
  
  # 6. NegBin-Q, "exact" vs. approximate estimates
  df_nbin_Q_approx <- tibble(
    Date = rep(incidence_subset$Date[t_ends], 2),
    R = c(R_hat$nbin_Q, R_hat$nbin_Q_approx),
    # Wald CI
    lwr = c(
      R_hat$nbin_Q - qnorm(0.975) * R_hat_sd$nbin_Q,
      R_hat$nbin_Q_approx - qnorm(0.975) * R_hat_sd$nbin_Q_approx
    ),
    upr = c(
      R_hat$nbin_Q + qnorm(0.975) * R_hat_sd$nbin_Q,
      R_hat$nbin_Q_approx + qnorm(0.975) * R_hat_sd$nbin_Q_approx
    ),
    Model = factor(rep(
      c("NegBin-Q", "NegBin-Q approx."),
      each = length(models_nbin_Q)
    ))
  )
  p_nbin_Q_exact_vs_approx <- ggplot(
    df_nbin_Q_approx,
    aes(
      x = Date, y = R, ymin = lwr, ymax = upr, color = Model, fill = Model
    )
  ) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_ribbon(color = NA, alpha = 0.5) +
    scale_color_manual(
      name = "Model",
      values = c(model_colors["NegBin-Q"], "NegBin-Q approx." = "gray30"),
      labels = c("NegBin-Q" = "NegBin-Q\n", "NegBin-Q approx." = "NegBin-Q\napprox.")
    ) +
    scale_fill_manual(
      name = "Model",
      values = c(model_colors["NegBin-Q"], "NegBin-Q approx." = "gray30"),
      labels = c("NegBin-Q" = "NegBin-Q\n", "NegBin-Q approx." = "NegBin-Q\napprox.")
    ) +
    labs(
      title = "NegBin-Q vs. its approximation",
      y = expression(hat(R))
    ) +
    coord_cartesian(
      ylim = c(0, 3), 
      xlim = range(incidence_subset$Date)
    ) +
    plot_theme

  ret <- list(
    R_hat = R_hat, R_hat_sd = R_hat_sd, disp = disp, AIC = AIC_vals,
    plt = list(p_incidence = p_incidence, p_pois_vs_qpois = p_pois_vs_qpois, 
               p_nbin_L_vs_nbin_Q = p_nbin_L_vs_nbin_Q, p_disp = p_disp,
               p_nbin_L_vs_qpois = p_nbin_L_vs_qpois,
               p_nbin_Q_exact_vs_approx = p_nbin_Q_exact_vs_approx)
  )
  return(ret)
}
