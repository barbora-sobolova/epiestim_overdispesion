library("tidyverse")
library("patchwork")
library("EpiEstim")  # For the discr_si() function
library("gamlss")  # For the gamlss() function
theme_set(theme_bw())

# Script illustrating the contributions to the log-likelihood for the 
# Poisson model vs. Negative binomial models. The identity link is used in the
# estimation.

# Read the data ================================================================

# We use one estimation window from the influenza data
window_width <- 7

# Data are from:
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011439#pcbi.1011439.s001
incidence <- read_csv("data/flu/daily_flu.csv")

# Select one time window for the illustration
start_date <- as.Date("2010-01-03")
end_date <- as.Date("2010-01-09")

# Subset the incidence
incidence_subset <- incidence %>%
  filter(Date >= start_date & Date <= end_date) %>% 
  mutate(Day = 1:window_width)

# Calculate the serial interval distribution ===================================

# Parameters of SI distribution same as Nash et al 2023 paper
mean_si <- 3.6
std_si <- 1.6
si_distr <- discr_si(seq(0, nrow(incidence) - 1), mean_si, std_si)

# Calculate Lambdas by weighting the past incidence ============================

# Calculate Lambda for all time points (adapted from EpiEstim)
lambda <- vector(mode = "numeric", length = nrow(incidence))
lambda[1] <- NA
for (t in seq(2, nrow(incidence))) {
  lambda[t] <- sum(
    si_distr[seq_len(t)] * incidence$Cases[seq(t, 1)])
}
lambda_subset <- lambda[incidence$Date >= start_date & 
                          incidence$Date <= end_date]

# Create the pairs of the mean values and the observations
pairs <- cbind(lambda_subset, incidence_subset$Cases)
colnames(pairs) <- c("lambda", "Cases")
model_mat <- as.data.frame(pairs)

# Fit the models ===============================================================

# Poisson model
mod_pois <- glm(
  data = model_mat,
  Cases ~ lambda - 1, 
  family = poisson(link = "identity")
)
# NegBin-Q model
mod_nbin_Q <- gamlss(
  data = model_mat,
  formula = Cases ~ lambda - 1,
  family = NBI(mu.link = "identity", sigma.link = "log"), 
  trace = FALSE
)
# NegBin-L model
mod_nbin_L <- gamlss(
  data = model_mat,
  formula = Cases ~ lambda - 1,
  family = NBII(mu.link = "identity", sigma.link = "log"), 
  trace = FALSE
)

# Extract the results ==========================================================

# Extract the coefficients
R_hat <- list(
  pois = mod_pois$coefficients,
  nbin_Q = mod_nbin_Q$mu.coefficients,
  nbin_L = mod_nbin_L$mu.coefficients
)
# Extract the dispersion parameters
disp_pars <- list(
  nbin_Q = exp(mod_nbin_Q$sigma.coefficients[1]),
  nbin_L = exp(mod_nbin_L$sigma.coefficients[1])
)
# Extract the value of the loglikelihood from the AIC values
ML <- list(
  pois = (2 - mod_pois$aic) / 2,
  nbin_Q = (4 - mod_nbin_Q$aic) / 2,
  nbin_L = (4 - mod_nbin_L$aic) / 2
)

# Prepare data frames for plotting =============================================

# Sequences for plotting
R_seq_lgt <- 1000
R_seq <- seq(0.1, 4.6, length = R_seq_lgt)
means <- c(R_seq %*% t(lambda_subset))
obs_long_vec <- rep(incidence_subset$Cases, each = R_seq_lgt)

# Create the data frames with the likelihood curves per observation

# Poisson
df_pois_llik <- tibble(
  R = rep(R_seq, times = window_width),
  llik = dpois(obs_long_vec, means, log = TRUE),
  Day = factor(rep(1:window_width, each = R_seq_lgt)),
  log_likelihood_of = "Individual obs.",
  model = "Poisson"
) %>% group_by(Day) %>% 
  mutate(
    # Find the maximum value of the log-likelihood for each curve to shift its
    # peak to 0.
    max_llik = max(llik),
    rel_llik = llik - max_llik
  ) %>% ungroup()
df_pois_llik_total <- df_pois_llik %>% group_by(R) %>% 
  summarise(
    llik = sum(llik),
    max_llik = ML$pois,
    rel_llik = llik - max_llik,
    log_likelihood_of = "Total",
    Day = "Total",
    model = "Poisson"
  )
# NegBin-Q
df_nbin_Q_llik <- tibble(
  R = rep(R_seq, times = window_width),
  llik = dnbinom(
    obs_long_vec, 
    mu = means, 
    size = 1 / disp_pars$nbin_Q,  # Same for all curves
    log = TRUE
  ),
  Day = factor(rep(1:window_width, each = R_seq_lgt)),
  log_likelihood_of = "Individual obs.",
  model = "NegBin-Q"
) %>% group_by(Day) %>% 
  mutate(
    # Find the maximum value of the log-likelihood for each curve to shift its
    # peak to 0.
    max_llik = max(llik),
    rel_llik = llik - max_llik
  ) %>% ungroup()
df_nbin_Q_llik_total <- df_nbin_Q_llik %>% group_by(R) %>% 
  summarise(
    llik = sum(llik),
    max_llik = ML$nbin_Q,
    rel_llik = llik - max_llik,
    log_likelihood_of = "Total",
    Day = "Total",
    model = "NegBin-Q"
  )
# NegBin-L
df_nbin_L_llik <- tibble(
  R = rep(R_seq, times = window_width),
  llik = dnbinom(
    obs_long_vec, 
    mu = means, 
    size = means / disp_pars$nbin_L,  # Dispersion parameter same for all curves
    log = TRUE
  ),
  Day = factor(rep(1:window_width, each = R_seq_lgt)),
  log_likelihood_of = "Individual obs.",
  model = "NegBin-L"
) %>% group_by(Day) %>% 
  mutate(
    # Find the maximum value of the log-likelihood for each curve to shift its
    # peak to 0.
    max_llik = max(llik),
    rel_llik = llik - max_llik
  ) %>% ungroup()
df_nbin_L_llik_total <- df_nbin_L_llik %>% group_by(R) %>% 
  summarise(
    llik = sum(llik),
    max_llik = ML$nbin_L,
    rel_llik = llik - max_llik,
    log_likelihood_of = "Total",
    Day = "Total",
    model = "NegBin-L"
  )

# Prepare the data frame with the point estimates
df_R_hat <- data.frame(
  R_hat = unlist(R_hat),
  model = factor(
    c("Poisson", "NegBin-Q", "NegBin-L"),
    levels = c("Poisson", "NegBin-L", "NegBin-Q")
  ),
  vline = "R_ML"
)
df_R_1 <- df_R_hat |> mutate(
  R_hat = 1,
  vline = "R_1"
)
df_plot_vlines <- rbind(df_R_hat, df_R_1)

df_llik <- rbind(
  df_pois_llik, df_nbin_L_llik, df_nbin_Q_llik,
  df_pois_llik_total, df_nbin_L_llik_total, df_nbin_Q_llik_total
) %>% 
  mutate(model = factor(model, levels = c("Poisson", "NegBin-L", "NegBin-Q")))

# Generate the black-and-white illustrative plot ===============================

# The incidence plot with Lambda_t
p_incidence_bw <- ggplot() +
  geom_segment(
    data = incidence_subset, 
    aes(x = Date, xend = Date, y = 0, yend = Cases)
  ) +
  geom_vline(aes(xintercept = start_date, linetype = "Cases"), alpha = 0) +
  geom_line(
    aes(x = incidence_subset$Date, y = lambda_subset, color = "lambda"),
    key_glyph = "timeseries"
  ) +
  scale_color_manual(
    values = "gray", 
    labels = expression(Lambda[t]),
    name = ""
  ) +
  scale_linetype_manual(
    values = "solid", 
    name = "",
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  guides(color = guide_legend(ncol = 2)) +
  labs(y = "Cases", color = "Day") +
  scale_x_date(date_labels = "%d %b ")

# 3 facets of the individual likelihood curves
p_loglik_individual_bw <- ggplot() +
  geom_line(
    data = df_llik, 
    aes(
      x = R, 
      y = rel_llik, 
      group = Day, 
      alpha = log_likelihood_of, 
      linewidth = log_likelihood_of
    )
  ) +
  geom_vline(
    data = df_plot_vlines, 
    aes(xintercept = R_hat, linetype = vline),
    linewidth = 0.3
  ) +
  coord_cartesian(ylim = c(-15, 0)) +
  scale_linetype_manual(
    values = c("R_ML" = "dashed", "R_1" = "dotted"),
    labels = c("R_ML" = expression(hat(R)[ML]), "R_1" = "R = 1"),
    name = ""
  ) +
  scale_linewidth_manual(
    values = c("Total" = 0.4, "Individual obs." = 0.3), 
    guide = "none"
  ) +
  scale_alpha_manual(
    values = c("Individual obs." = 0.3, "Total" = 1), 
    name = "Relative\nlog-likelihood",
    labels = c("Individual \nobs.", "Total")
  ) +
  labs(
    y = "Relative log-likelihood",
    x = "R"
  ) +
  facet_wrap(~model, nrow = 3) +
  theme(strip.background = element_blank(), legend.key.spacing.y = unit(0.3, "cm"))
p_incidence_loglik_bw <- p_incidence_bw + p_loglik_individual_bw + 
  plot_layout(design = "A\nB\nB\nB\nB")

ggsave("figure/Incidence_and_llik_bw_identity.pdf", p_incidence_loglik_bw, width = 6, 
       height = 5)
