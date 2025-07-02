library(tidyverse)
library(ggplot2)
library(EpiEstim)
library(gamlss)
library(patchwork)
source("scripts/analyse_Rt.R")

# Running the analysis on 3 datasets:
# - flu among USA military personnel
# - COVID-19 in Austria
# - Ebola in Guinea

# Set the locale to English
Sys.setlocale("LC_ALL","English")

# List for storing the parameters for the 3 datasets
params <- replicate(3, list())
names(params) <- c("flu", "covid", "ebola")

# Read the flu data --------------------------------------------------

# Data are from:
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011439#pcbi.1011439.s001
incidence_flu <- read_csv("data/flu/daily_flu.csv")

# NOTES:
# Window size used is 7 days
# SI distribution is from the paper below (mean 3.6 days, sd 1.6 days)

# We run the analysis on the whole set
params$flu$incidence <- incidence_flu
params$flu$start_date <- incidence_flu$Date[1]
params$flu$end_date <- incidence_flu$Date[nrow(incidence_flu)]
params$flu$window_width <- 7
# Parameters of SI distribution same as Nash et al 2023 paper
params$flu$mean_si <- 3.6
params$flu$std_si <- 1.6

# Read the covid data --------------------------------------------------

# Data from the ECDC web:
# https://www.ecdc.europa.eu/en/publications-data/data-daily-new-cases-covid-19-eueea-country
incidence_covid <- read_csv("data/covid/covid_ecdc.csv") %>%
  filter(geoId == "AT") %>%
  dplyr::select(dateRep, cases) %>%
  rename("Date" = "dateRep", "Cases" = "cases") %>%
  mutate(
    Date = as.Date(Date, format = "%d/%m/%Y")
  )
# Order the data chronologically
incidence_covid <- incidence_covid[nrow(incidence_covid):1, ]

# We run the analysis on a subset just to keep it shorter
params$covid$incidence <- incidence_covid
params$covid$start_date <- as.Date("2021-10-01")
params$covid$end_date <- as.Date("2022-03-31")
params$covid$window_width <- 13
# Parameters of the serial interval distribution from (Richter et al, 2022)
params$covid$mean_si <- 3.37
params$covid$std_si <- 1.83

# Read the ebola data --------------------------------------------------

# Data are from:
# https://github.com/fintzij/stemr/tree/master/data
incidence_ebola_daily <-  read_csv(file = "data/ebola/green2022_data_frame.csv")

# Rescale the parameters to the weekly incidence
# Parameters of SI distribution taken from Green et al 2022
mu <- 15.3
sd <- 9.3
# shape
alpha <- (mu / sd)^2
# rate
beta <- (mu / sd^2)
# divide by 7 means
scaled_beta <- beta / (1 / 7)
scaled_mean <- alpha / scaled_beta
scaled_sd <- sqrt(alpha / (scaled_beta^2))

# Aggregate the data into weekly format
incidence_ebola_weekly <- incidence_ebola_daily %>%
  group_by(week = lubridate::floor_date(date, "week")) %>%
  summarise(Cases = sum(count)) %>%
  ungroup() %>%
  mutate(week = as.Date(week)) %>%
  filter(week > "2014-02-23") %>%
  mutate(Date = week)

params$ebola$incidence <- incidence_ebola_weekly
params$ebola$start_date <- as.Date("2014-03-30")
params$ebola$end_date <- max(incidence_ebola_weekly$Date)
params$ebola$window_width <- 4
params$ebola$mean_si <- scaled_mean
params$ebola$std_si <- scaled_sd

# Run the models --------------------------------------------------

results <- vector("list", 3)
names(results) <- names(params)
for (disease in names(results)) {
  results[[disease]] <- with(
    params[[disease]], 
    analyse_Rt(incidence, start_date, end_date, window_width, mean_si, std_si)
  )
}

# Plot the results --------------------------------------------------

# For Ebola, we have trouble fitting the NegBin-L model for 
# 19. and 26. April 2015. In these 2 estimation windows, the counts are 
# underdispersed. This can be seen also in the Poisson vs. quasi-poisson plot,
# where the dispersion parameter is estimated as < 1 and the resulting 
# confidence intervals are actually higher for Poisson.
# As a result, NegBin-L is trying to estimate its dispersion parameter as very 
# low and runs into the problem of numerical instabilities. For this reason we 
# replace these 2 NegBin-L estimates by the Poisson estimates.
# The changes are carried out in the plot only! If we want to report other 
# quantities, we would have to change them too.
dates_underdisp <- as.Date(c("2015-04-26", "2015-04-19"))
df_replace <- results$ebola$plt$p_nbin_L_vs_nbin_Q$data |> filter(
  Date %in% dates_underdisp & Model == "Poiss"
) |> mutate(Model = "NegBin-L")
df_updated <- rows_update(
  results$ebola$plt$p_nbin_L_vs_nbin_Q$data, 
  df_replace, 
  by = c("Date", "Model")
)
results$ebola$plt$p_nbin_L_vs_nbin_Q <- results$ebola$plt$p_nbin_L_vs_nbin_Q %+% df_updated
df_replace <- results$ebola$plt$p_nbin_L_vs_qpois$data |> filter(
  Date %in% dates_underdisp & Model == "Poiss"
) |> mutate(Model = "NegBin-L")
df_updated <- rows_update(
  results$ebola$plt$p_nbin_L_vs_qpois$data, 
  df_replace, 
  by = c("Date", "Model")
)
results$ebola$plt$p_nbin_L_vs_qpois <- results$ebola$plt$p_nbin_L_vs_qpois %+% df_updated

# Adjust individual plots
plots_to_adjust_date <- c("p_incidence", "p_nbin_L_vs_nbin_Q", "p_pois_vs_qpois",
                          "p_disp")
for (plt in plots_to_adjust_date) {
  results$flu$plt[[plt]] <- results$flu$plt[[plt]] +
    scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d")
}
for (plt in plots_to_adjust_date) {
  results$covid$plt[[plt]] <- results$covid$plt[[plt]] +
    scale_x_date(date_breaks = "3 months", date_labels = "%b %Y")
}
for (plt in plots_to_adjust_date) {
  results$ebola$plt[[plt]] <- results$ebola$plt[[plt]] +
    scale_x_date(date_breaks = "4 months", date_labels = "%b %Y")
}
  
# Prepare the annotations
annotations <- list()
annotations$flu <- plot_annotation(
  title = "Influenza,\nUSA Active Military Personnel,\n2009-2010",
  theme = theme(
    plot.title = element_text(
      size = 16,
      hjust = 0.5,
      face = "bold",
      lineheight = 0.9
    )
  )
)
annotations$covid <- plot_annotation(
  title = "COVID-19,\nAustria,\n2021-2022",
  theme = theme(
    plot.title = element_text(
      size = 16,
      hjust = 0.5,
      face = "bold",
      lineheight = 0.9
    )
  )
)
annotations$ebola <- plot_annotation(
  title = "Ebola,\nGuinea,\n2014-2015",
  theme = theme(
    plot.title = element_text(
      size = 16,
      hjust = 0.5,
      face = "bold",
      lineheight = 0.9
    )
  )
)

# Compose the plots
p_legend <- wrap_elements(ggpubr::get_legend(results$flu$plt$p_nbin_L_vs_nbin_Q))
composite_plot_elements <- overdisp_plot_elements <- vector("list", 3)
names(composite_plot_elements) <- names(overdisp_plot_elements) <- names(params)
for (disease in names(composite_plot_elements)) {
  composite_plot_elements[[disease]] <- (
    with(
      results[[disease]]$plt,
      p_incidence / p_pois_vs_qpois / p_nbin_L_vs_nbin_Q
    ) +
      plot_layout(heights = c(1, 1, 1), guides = "collect") &
      theme(legend.position = "none") & 
      annotations[[disease]]
  ) |>
    wrap_elements()
  
  overdisp_plot_elements[[disease]] <- (
    with(
      results[[disease]]$plt,
      p_incidence / p_disp
    ) +
      plot_layout(heights = c(1, 1), guides = "collect") &
      theme(legend.position = "none") & 
      annotations[[disease]]
  ) |>
    wrap_elements()
}

composite_plot <- (
  with(
    composite_plot_elements, 
    flu | covid | ebola | p_legend
  )
) +
  plot_layout(widths = c(3, 3, 3, 1))
ggsave("figure/composite_plot_overdisp_longer.pdf", composite_plot, width = 14, height = 11)

overdisp_plot <- (
  with(
    overdisp_plot_elements, 
    flu | covid | ebola | p_legend
  )
)
ggsave("figure/overdisp_longer_plot.pdf", overdisp_plot, width = 14, height = 6.5)

# Supplementary plots --------------------------------------------------

# NegBin-L vs. quasi-Poisson
p_nbin_L_legend <- wrap_elements(ggpubr::get_legend(results$flu$plt$p_nbin_L_vs_qpois))
nbin_L_plot_elements <- vector("list", 3)
names(nbin_L_plot_elements) <- names(params)
for (disease in names(nbin_L_plot_elements)) {
  nbin_L_plot_elements[[disease]] <- (
    with(
      results[[disease]]$plt,
      p_incidence / p_nbin_L_vs_qpois
    ) +
      plot_layout(heights = c(1, 1), guides = "collect") &
      theme(legend.position = "none") & 
      annotations[[disease]]
  ) |>
    wrap_elements()
}

nbin_L_plot <- (
  with(
    nbin_L_plot_elements, 
    flu | covid | ebola | p_nbin_L_legend
  )
) +
  plot_layout(widths = c(3, 3, 3, 1))
ggsave("figure/nbin_L_vs_qpis_plot.pdf", nbin_L_plot, width = 14, height = 6.5)

# NegBin-Q "exact" vs. approximate
p_nbin_Q_legend <- wrap_elements(ggpubr::get_legend(results$flu$plt$p_nbin_Q_exact_vs_approx))
nbin_Q_plot_elements <- vector("list", 3)
names(nbin_Q_plot_elements) <- names(params)
for (disease in names(nbin_Q_plot_elements)) {
  nbin_Q_plot_elements[[disease]] <- (
    with(
      results[[disease]]$plt,
      p_incidence / p_nbin_Q_exact_vs_approx
    ) +
      plot_layout(heights = c(1, 1), guides = "collect") &
      theme(legend.position = "none") & 
      annotations[[disease]]
  ) |>
    wrap_elements()
}

nbin_Q_plot <- (
  with(
    nbin_Q_plot_elements, 
    flu | covid | ebola | p_nbin_Q_legend
  )
) +
  plot_layout(widths = c(3, 3, 3, 1))
ggsave("figure/nbin_Q_approximation_plot.pdf", nbin_Q_plot, width = 14, height = 6.5)
