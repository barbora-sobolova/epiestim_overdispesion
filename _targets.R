library(targets)
library(tarchetypes)
library(qs2)

# Set targets options
tar_option_set(
  packages = c(
    "ggplot2",
    "patchwork",
    "purrr",
    "here",
    "tidyr",
    "dplyr",
    "qs2",
    "EpiEstim",
    "gamlss",
    "readr"
  ),
  format = "qs", # Use qs format (qs2 is used via repository option)
  memory = "transient", # Free memory after each target completes
  garbage_collection = TRUE, # Run garbage collection
  repository = "local", # Use qs2 backend for storage
  error = "continue" # Continue pipeline when targets fail
)
tar_source(files = "functions")

# set the ggplot theme
ggplot2::theme_set(ggplot2::theme_bw())

# Set the global objects =======================================================
model_colors <- c(
  "Poiss" = "#009E73",
  "Q-Poiss" = "#F0E442",
  "NegBin-L" = "#56B4E9",
  "NegBin-Q" = "#CC79A7"
)

# Parameters common for all simulation runs
global_params <- list(
  mean_si = 7.5,
  std_si = 2.1,
  n_init = 14,
  short_window = 7,
  long_window = 14,
  n_sim = 1000L,
  # We don't use the weekday effects in the simulation. However, they are left
  # here just in case we need them in the future.
  weekday = c(1.05, 1.40, 1.75, 1.40, 1.05, 0.21, 0.14),
  base_seed = 9786L
)

# Size of the final plot
plot_size <- list(width = 11, height = 11.7)

# If we plot only half of the scenarios, by how much do we divide the height?
plot_halving_coeff <- 1.75

# Scenario blocks for the static branching:
# NegBin-L
# NegBin-Q
# Poisson
# None of them involves weekday effects, but they can be added by using
# "weekday_yes"
outer_scenarios <- data.frame(
  distribution = c("NegBin-L", "NegBin-Q", "Poiss"),
  weekday_effect = c("weekday_no", "weekday_no", "weekday_no")
) |>
  dplyr::mutate(
    scenario_id = paste(distribution, weekday_effect, sep = "_")
  )

# Define pipeline ==============================================================
list(
  # Serial interval, common for all scenarios
  tar_target(
    si,
    with(global_params, discr_si(seq_len(n_init), mean_si, std_si))
  ),
  # Static branching over 4 scenario blocks defined in `outer_scenarios`
  tar_map(
    unlist = FALSE,
    values = outer_scenarios,
    names = scenario_id,
    # Get the data frame with scenario parameters
    tar_target(scenarios, create_scenario_grid(distribution)),
    # Set the weekday effect. If the weekday effect is not present, we multiply
    # the simulated mean value by 1 under the hood, so there is indeed no
    # effect.
    tar_target(
      weekday_effect_vector,
      if (weekday_effect == "weekday_yes") {
        global_params$weekday
      } else {
        rep(1, global_params$n_init)
      }
    ),
    # Initial values. Sample iid counts using a reasonable data generating
    # mechanism.
    tar_target(
      init,
      initialize_trajectory(
        global_params$n_init,
        scenarios$init_magnitude,
        scenarios$init_sd,
        weekday_effect = weekday_effect_vector,
        seed = global_params$base_seed + scenarios$init_seed
      ),
      pattern = map(scenarios),
      iteration = "list"
    ),
    # Simulations
    tar_target(
      trajectories,
      with(
        global_params,
        generate_trajectories(
          n_sim = n_sim,
          init,
          R_eff = scenarios$R_eff,
          si = si,
          lgt = n_init + long_window,
          model = distribution,
          nb_size = scenarios$nb_size,
          seed = global_params$base_seed + scenarios$scenario_number,
          weekday_effect = weekday_effect_vector
        )
      ),
      pattern = map(init, scenarios),
      iteration = "list"
    ),
    # Estimation
    tar_target(
      df_R_hat_raw,
      create_results_df(
        trajectories$X,
        trajectories$Lambda,
        global_params$short_window,
        global_params$long_window
      ),
      pattern = map(trajectories),
      iteration = "list"
    ),
    # Replace the divergent NegBin estimates by the Poisson ones
    tar_target(
      df_R_hat,
      replace_divergent(df_R_hat_raw),
      pattern = map(df_R_hat_raw),
      iteration = "list"
    ),
    # Plot the trajectories and the coverage
    tar_target(
      plot_panels,
      list(
        trajectories = plot_trajectories(
          trajectories$X,
          global_params$short_window,
          global_params$n_init
        ),
        coverage = plot_coverage(
          scenarios$R_eff,
          scenarios$nb_size,
          df_R_hat,
          seq(0, 1, by = 0.01),
          distribution,
          weekday_effect,
          model_colors
        ),
        meta = plot_metadata(
          scenarios$R_eff,
          scenarios$nb_size,
          scenarios$magnitude,
          distribution
        )
      ),
      pattern = map(trajectories, df_R_hat, scenarios),
      iteration = "list"
    ),
    # How many did converge?
    tar_target(
      convergence_tables,
      summarize_convergence(
        df_R_hat,
        scenarios$magnitude,
        scenarios$nb_size,
        scenarios$R_eff,
        # Model name from the outer scenarios data frame
        distribution
      ),
      pattern = map(df_R_hat, scenarios),
      iteration = "list"
    ),
    # Save the convergence table
    tar_target(
      saved_tables,
      {
        convergence <- convergence_tables |>
          purrr::map("df_convergence") |>
          dplyr::bind_rows()
        unstable <- convergence_tables |>
          purrr::map("df_unstable") |>
          dplyr::bind_rows()
        write_csv(
          convergence,
          paste0("tables/", scenario_id, "_convergence.csv")
        )
        write_csv(unstable, paste0("tables/", scenario_id, "_unstable.csv"))
      }
    ),
    # Plot the density of R estimates and its SEs
    tar_target(
      plot_density_panels,
      plot_dens(
        df_R_hat,
        scenarios$R_eff,
        1 / scenarios$nb_size,
        scenarios$magnitude,
        distribution,
        model_colors
      ),
      pattern = map(df_R_hat, scenarios),
      iteration = "list"
    ),
    # Save plots
    tar_target(saved_figures, {
      if (distribution == "Poiss") {
        # For Poisson, we have only half the scenarios as for the rest, so the
        # height of the resulting plot must be divided by 2. We divide by less
        # than 2 to allow for some space for the title.
        plot_height <- plot_size$height / plot_halving_coeff
      } else {
        plot_height <- plot_size$height
      }

      # Save the coverage plot
      p_simulation <- compose_coverage_patches(
        plot_panels,
        global_params$short_window,
        global_params$long_window
      )
      save_plot(
        p_simulation,
        paste(scenario_id, "simulation_coverage", sep = "_"),
        width = plot_size$width,
        height = plot_height
      )
      # Save the distribution of the R estimates
      p_densities <- compose_dens_patches(
        # Extract only the R estimates
        list(
          R_hat = purrr::map(plot_density_panels, "R_hat"),
          se_hat = purrr::map(plot_density_panels, "se_hat")
        ),
        purrr::map(plot_panels, "meta"),
        global_params$short_window,
        global_params$long_window
      )
      save_plot(
        p_densities,
        paste(scenario_id, "Rhat_density", sep = "_"),
        width = plot_size$width,
        height = plot_height
      )
    }),
    tar_target(
      saved_figures_poster,
      {
        # Create the coverage plot for a subset of scenarios. We will use only
        # NegBin-L and NegBin-Q versions.
        # CURRENTLY, THIS PLOT NEEDS ADJUSTMENTS REGARDING THE TEXT SIZES.
        # NEEDS TO BE RESOLVED BEFORE GENERATING THE PLOTS FOR THE POSTER.
        generate_poster_figure <- distribution == "NegBin-Q" ||
          (distribution == "NegBin-L" && weekday_effect == "weekday_no")
        if (generate_poster_figure) {
          p_simulation_poster <- compose_coverage_patches(
            plot_panels[scenarios$R_eff == 1.5 & scenarios$magnitude == "low"],
            global_params$short_window,
            global_params$long_window,
            panel_widths = c(1.3, 2, 3)
          ) &
            theme(plot.margin = unit(c(5.5, 5.5, 0.5, 5.5), "points"))
          save_plot(
            p_simulation_poster,
            paste(scenario_id, "simulation_coverage_poster", sep = "_"),
            # Figure dimensions must be set manually to fit the poster page
            width = plot_size$width / 1.4,
            height = (plot_size$height / 1.4) / 2.35
          )
        }
      }
    )
  ),
  tar_target(saved_overdisp_est_plots, {
    overdisp_panels <- list(
      NegBin.L_weekday_no = purrr::map(
        plot_density_panels_NegBin.L_weekday_no,
        "overdisp_hat"
      ),
      NegBin.Q_weekday_no = purrr::map(
        plot_density_panels_NegBin.Q_weekday_no,
        "overdisp_hat"
      )
    )
    meta_panels <- list(
      NegBin.L = purrr::map(plot_panels_NegBin.L_weekday_no, "meta"),
      NegBin.Q = purrr::map(plot_panels_NegBin.Q_weekday_no, "meta")
    )
    p_overdisp <- compose_overdisp_patches(overdisp_panels, meta_panels)
    save_plot(
      p_overdisp,
      "overdisp_estimates",
      width = plot_size$width,
      height = plot_size$height
    )
  })
)
