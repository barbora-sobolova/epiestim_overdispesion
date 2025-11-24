#' Returns the theme of the composite plot legend
#'
#' @description This macro defines the font sizes used in the legend of the
#'   composite plots.
#' @return a ggplot theme
get_legend_theme <- function() {
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
}

#' Returns the x-axis limits for the R estimate distribution plot
#'
#' @description This function returns the x-axis limits for the plot of R
#'   estimate distribution, the plot of its standard errors and the plot of the
#'   overdispersion parameter estimate density. The limits are determined
#'   manually and returned as a list.
#' @param R_eff positive real, the true value of the effective reproductive
#'   number. The x-axis is centered around this value.
#' @param overdisp_true the true value of the overdispersion parameter used to
#'   generate the trajectory
#' @param magnitude string indicating whether the scenario operates with low, or
#'   high magnitudes.
#' @param distribution the true underlying count distribution
#' @return a list of 3 vectors of 2 elements - the limits for R estimate,
#'   the limits for its standard errors and the limits of the overdispersion
#'   estimates
get_xlim <- function(R_true, overdisp_true, magnitude, distribution) {
  if (distribution == "NegBin-L") {
    se_limits <- if (magnitude == "high") {
      c(0, 0.15)
    } else if (magnitude == "low") {
      c(0, 0.6)
    }
    overdisp_limits <- c(-1, min(10, 5 * overdisp_true))
  } else if (distribution == "NegBin-Q") {
    se_limits <- if (magnitude == "high") {
      c(0, 0.2)
    } else if (magnitude == "low") {
      c(0, 0.4)
    }
    overdisp_limits <- c(0, overdisp_true + 0.08)
  } else if (distribution == "Poiss") {
    se_limits <- if (magnitude == "high") {
      c(0, 0.05)
    } else if (magnitude == "low") {
      c(0, 0.3)
    }
    overdisp_limits <- c(NA, NA)
  }
  list(
    R_hat = R_true + c(-0.8, 0.8),
    se_hat = se_limits,
    overdisp_hat = overdisp_limits
  )
}

#' Plots density distributions of the estimates
#'
#' @description This function plots density distributions of the R estimates,
#'   highlighting the true value, and density distributions of the standard
#'   error estimates
#' @param df_R_hat a data frame with R estimates, after replacing the divergent
#'   runs by NAs, containing columns
#' @param R_true the true value of the effective reproduction number used to
#'   generate the trajectory
#' @param overdisp_true the true value of the overdispersion parameter used to
#'   generate the trajectory
#' @param magnitude string indicating whether the scenario operates with low, or
#'   high magnitudes.
#' @param dist_true string, the true count distribution of the data generating
#'   process
#' @param model_colors a named vector specifying the model colors. Must have
#'   names "Poiss", "Q-Poiss", "NegBin-L" and "NegBin-Q"
#' @return a nested list containing density plots: outer list has elements
#'   \code{R_hat} and \code{se_hat}, each containing a list of two plots
#'   (one per window length), and \code{overdisp_hat} containing just one
#'   density plot
plot_dens <- function(
  df_R_hat,
  R_true,
  overdisp_true,
  magnitude,
  dist_true,
  model_colors
) {
  # Create a denser grid where we want `geom_line(stat = "density")` to
  # calculate the density estimate. the default 512 occasionally creates a
  # ragged line.
  density_estimation_grid <- 1024

  # Split the data frames to create 2 separate plots for each window width. This
  # way, it's easier to put the final plot together from multiple blocks
  df_R_hat_split <- split.data.frame(
    df_R_hat,
    df_R_hat$window_len_fct
  )

  # Calculate the x-limits for the plot of R estimates, its standard errors and
  # the dispersion parameter estimates
  limits_x <- get_xlim(R_true, overdisp_true, magnitude, dist_true)

  p_R_hat <- p_se_hat <- p_overdisp_hat <- vector("list", 2)
  names(p_R_hat) <- names(p_se_hat) <- names(df_R_hat_split)
  for (k in 1:2) {
    # Plot the point estimates of R_eff (R_hat)
    p_R_hat[[k]] <- ggplot(df_R_hat_split[[k]], aes(x = R, color = model)) +
      geom_line(
        stat = "density",
        linewidth = 1,
        alpha = 0.6,
        na.rm = TRUE,
        n = density_estimation_grid
      ) +
      geom_vline(aes(xintercept = R_true, linetype = "R_true"), color = "red") +
      scale_color_manual(values = model_colors) +
      scale_linetype_manual(
        values = c("R_true" = "dashed"),
        labels = c("R_true" = expression(R[t]))
      ) +
      labs(
        x = expression(hat(R)),
        color = "Model",
        linetype = NULL,
        y = "density"
      ) +
      # A couple of values might get clipped.
      coord_cartesian(xlim = limits_x$R_hat)

    # Plot the standard errors of the R_eff estimates
    p_se_hat[[k]] <- ggplot(df_R_hat_split[[k]], aes(x = se, color = model)) +
      geom_line(
        stat = "density",
        linewidth = 1,
        alpha = 0.6,
        na.rm = TRUE,
        bounds = c(0, Inf),
        # Increase the grid density a little more for the Poison SEs, where
        # the line is more rough than for the rest.
        n = if (dist_true == "Poiss") {
          2 * density_estimation_grid
        } else {
          density_estimation_grid
        }
      ) +
      scale_color_manual(values = model_colors) +
      labs(
        x = expression(widehat(se)(hat(R))),
        color = "Model",
        y = "density"
      ) +
      # A couple of values might get clipped.
      coord_cartesian(xlim = limits_x$se_hat)
  }

  # Plot the point estimates of the overdispersion parameter. We use 2 facets
  # for 2 window lengths, therefore, we create the plot outside the for loop.
  if (dist_true == "Poiss") {
    # For the Poisson ground truth we don't plot anything
    p_overdisp_hat <- NULL
  } else {
    # We plot only parameters on the same scale. If the ground truth is
    # NegBin-Q, then only the NegBin-Q parameter estimates are on the same
    # scale. For NegBin-L, we can plot both the NegBin-L and the quasi-Poisson
    # estimate.
    if (dist_true == "NegBin-Q") {
      # set the x-axis label to psi for NegBin-Q
      x_label <- expression(psi)
      # Plot the geometries
      p_overdisp_hat <- ggplot() +
        geom_line(
          data = filter(df_R_hat, model == "NegBin-Q"),
          aes(x = overdisp, color = "NegBin-Q", linetype = window_len_fct),
          stat = "density",
          linewidth = 1,
          alpha = 0.6,
          na.rm = TRUE,
          bounds = c(0, Inf),
          n = density_estimation_grid
        )
    } else if (dist_true == "NegBin-L") {
      # set the x-axis label to xi for NegBin-L
      x_label <- expression(xi)
      # Plot the geometries
      p_overdisp_hat <- ggplot() +
        geom_line(
          data = filter(df_R_hat, model == "Q-Poiss"),
          # Shift the overdispersion parameter by 1 to make it comparable with
          # NegBin-L
          aes(x = overdisp - 1, color = "Q-Poiss", linetype = window_len_fct),
          stat = "density",
          linewidth = 1,
          alpha = 0.6,
          na.rm = TRUE,
          # Where the density shall be calculated. We need to plot quasi-Poisson
          # and NegBin-L separately to be able to set different bounds.
          bounds = c(-1, Inf),
          n = density_estimation_grid
        ) +
        geom_line(
          data = filter(df_R_hat, model == "NegBin-L"),
          aes(x = overdisp, color = "NegBin-L", linetype = window_len_fct),
          stat = "density",
          linewidth = 1,
          alpha = 0.6,
          na.rm = TRUE,
          bounds = c(0, Inf),
          n = density_estimation_grid
        ) +
        # A dummy line to achieve including the "NegBin-Q" label in the legend
        geom_line(
          aes(x = c(-10, -11), y = c(0, 0), color = "NegBin-Q"),
          linewidth = 1,
          alpha = 0.6
        )
    }

    # Extract the lengths of the estimation window to include it in plot
    # labels
    window_lengths <- unique(df_R_hat$window_len)

    p_overdisp_hat <- p_overdisp_hat +
      geom_vline(
        aes(xintercept = overdisp_true, linetype = "overdisp_true"),
        color = "red"
      ) +
      scale_color_manual(
        values = model_colors,
        breaks = c("Q-Poiss", "NegBin-L", "NegBin-Q")
      ) +
      scale_linetype_manual(
        name = NULL,
        values = c(
          "overdisp_true" = "dotted",
          "short" = "dotted",
          "long" = "solid"
        ),
        labels = c(
          "overdisp_true" = "overdispersion\nparameter",
          "short" = paste0(min(window_lengths), "-day window"),
          "long" = paste0(max(window_lengths), "-day window")
        )
      ) +
      labs(
        x = x_label,
        color = "Model",
        y = "density"
      ) +
      # A couple of values might get clipped.
      coord_cartesian(xlim = limits_x$overdisp_hat)
  }

  # Return the list
  list(R_hat = p_R_hat, se_hat = p_se_hat, overdisp_hat = p_overdisp_hat)
}

#' Plots nominal vs empirical coverage
#'
#' @description This function plots nominal vs. empirical coverage levels with
#'   a 45-degree reference line, faceted by model.
#' @param R_eff the true value of the effective reproduction number used to
#'   generate the trajectory
#' @param nb_size positive real value, the size parameter of the negative
#'   binomial distribution for "NegBin-Q" and "NegBin-L"
#' @param df_R_hat a data frame with R estimates, after replacing the divergent
#'   runs by NAs, containing columns
#'   \code{R_hat}, \code{se_hat} and \code{converged}
#' @param window the length of the estimation window, can be dropped when we get
#'   rid of the normal approximation part
#' @param nominal_covr a vector of the nominal coverage levels
#' @return a list containing 2 ggplots of coverage
plot_coverage <- function(
  R_eff,
  nb_size,
  df_R_hat,
  nominal_covr,
  distribution,
  weekday_effect,
  model_colors
) {
  df_coverage <- create_coverage_df(
    R_eff,
    nb_size,
    df_R_hat,
    nominal_covr,
    distribution,
    weekday_effect
  )

  # Split the data frames to create 2 separate plots instead of facets. This
  # way, it's easier to put the final plot together from multiple blocks
  df_split <- split.data.frame(
    df_coverage,
    df_coverage$window_len_fct
  )

  # Don't show the legend for the theoretical Poisson coverage for the
  # "NegBin-Q" distribution
  if (distribution == "NegBin-L") {
    linetype_guide <- "legend"
  } else {
    linetype_guide <- "none"
  }

  p_coverage <- vector("list", 2)
  names(p_coverage) <- names(df_split)
  for (k in 1:2) {
    p_coverage[[k]] <- ggplot() +
      geom_abline(intercept = 0, slope = 1, color = "grey40") +
      geom_line(
        data = df_split[[k]],
        mapping = aes(
          x = covr_nominal,
          y = covr_empirical,
          color = model,
          linetype = type,
          linewidth = type,
          alpha = type
        )
      ) +
      scale_color_manual(
        values = c(model_colors, "theoretical Poisson" = "black"),
        breaks = names(model_colors),
        labels = names(model_colors)
      ) +
      scale_linetype_manual(
        values = c(
          "Empirical coverage" = "solid",
          "theoretical Poisson coverage" = "dashed"
        ),
        breaks = "theoretical Poisson coverage",
        labels = c(
          "theoretical Poisson coverage" = "theoretical\nPoisson\ncoverage"
        ),
        guide = linetype_guide
      ) +
      scale_linewidth_manual(
        values = c(
          "Empirical coverage" = 1.5,
          "theoretical Poisson coverage" = 0.5
        ),
        guide = "none"
      ) +
      scale_alpha_manual(
        values = c(
          "Empirical coverage" = 0.6,
          "theoretical Poisson coverage" = 1
        ),
        guide = "none"
      ) +
      scale_x_continuous(
        breaks = seq(0, 1, by = 0.2),
        minor_breaks = seq(0, 1, by = 0.1)
      ) +
      scale_y_continuous(
        breaks = seq(0, 1, by = 0.2),
        minor_breaks = seq(0, 1, by = 0.1)
      ) +
      labs(
        x = "Nominal coverage",
        y = "Empirical coverage",
        color = "Model",
        shape = NULL,
        linetype = NULL
      )
  }
  p_coverage
}

#' Plots the fan-out trajectories
#'
#' @description This function plots the fan-out of the trajectory from its fixed
#'   initial values. The next step is plotting the trajectory for the longer
#'   estimation windows and making the part corresponding to the shorter window
#'   lighter/darker.
#' @param X a matrix of the simulated counts, one column per simulation run
#' @param short_window an integer, the length of the short estimation window
#' @param n_init an integer, the number of initial values of the trajectory
#' @return a ggplot object
plot_trajectories <- function(X, short_window, n_init) {
  df_trajectories <- reshape2::melt(
    # Display only the first 100 simulation runs to make the plot readable
    X[, seq_len(min(100, ncol(X)))],
    varnames = c("day", "trajectory"),
    value.name = "cases"
  ) |>
    mutate(
      segment = dplyr::case_when(
        day < n_init ~ "Initial",
        day >= n_init & day < short_window + n_init ~ "Short window",
        day >= short_window + n_init ~ "Long window"
      )
    )

  # Set the offset for the brackets
  max_cases <- max(df_trajectories$cases)
  bracket_offset <- max_cases * 0.08  # 8% of the max value

  # The long window is the length of the trajectory minus the initialization
  long_window <- nrow(X) - n_init

  ggplot() +
    geom_line(
      data = df_trajectories,
      mapping = aes(
        x = day,
        y = cases,
        group = trajectory,
        color = segment,
        alpha = segment
      )
    ) +
    ggpubr::geom_bracket(
      xmin = n_init + 1,
      xmax = nrow(X),
      y.position = max_cases + 2.5 * bracket_offset,
      label = paste0(long_window, "-day window"),
      label.size = 3
    ) +
    ggpubr::geom_bracket(
      xmin = n_init + 1,
      xmax = n_init + short_window,
      y.position = max_cases + 0.5 * bracket_offset,
      label = paste0(short_window, "-day window"),
      label.size = 3
    ) +
    scale_color_manual(
      values = c(
        "Initial" = "black",
        "Short window" = "black",
        "Long window" = "gray30"
      ),
      guide = "none"
    ) +
    scale_alpha_manual(
      values = c(
        "Initial" = 1,
        "Short window" = 0.1,
        "Long window" = 0.1
      ),
      guide = "none"
    ) +
    coord_cartesian(ylim = c(0, max_cases + 3.5 * bracket_offset)) +
    labs(x = "Day", y = "Cases") +
    theme(axis.title = element_text(size = 16))
}

#' Plots the metadata of the simulation scenario
#'
#' @description This function plots the information about the true parameter
#'   values used in the simulation scenario.
#' @param R_eff a positive real, the true value of the effective reproduction
#'   number
#' @param nb_size a positive real, the true value of the dispersion parameter
#' @param magnitude string values, either "high", or "low"
#' @param distribution string value indicating the count distribution:
#'   "NegBin-L", "NegBin-Q", or "Poiss"
#' @return a ggplot object
plot_metadata <- function(R_eff, nb_size, magnitude, distribution) {
  df_text <- data.frame(
    x = rep(1, 3),
    y = c(1, 0, -1),
    label = c(
      paste0("R[t] == ", R_eff),
      if (distribution == "NegBin-L") {
        paste0("xi == ", 1 / nb_size)
      } else if (distribution == "NegBin-Q") {
        paste0("psi == ", 1 / nb_size)
      } else {
        NA
      },
      paste0("Magn.: ", gsub("_.*", "", magnitude))
    )
  )
  # Remove the dispersion parameter, when it's not present for the Poisson
  # distribution.
  if (distribution == "Poiss") {
    df_text <- df_text |> filter(!is.na(label))
  }
  ggplot(df_text, aes(x = x, y = y, label = label)) +
    geom_text(hjust = 0, parse = TRUE) +
    coord_cartesian(ylim = c(-3, 3), xlim = c(0.997, 1.02)) +
    # Must specify the labels as an empty string, otherwise a later call of
    # theme() will override the labels as default "x" and "y"
    labs(x = "", y = "") +
    theme_void()
}

#' Compose coverage plot patches into the final plot
#'
#' @description This function composes coverage and trajectory plots from all
#'  the scenario runs together, including titles and labels, using the patchwork
#'  approach.
#' @param plot_panels a list of lists of patchwork patches. The number of
#'   simulation scenarios is the length of the outer list. Each inner list
#'   contains 4 plots - simulation trajectories, coverage for the short and
#'   long estimation window and a plot with parameter values.
#' @param short_window an integer, the length of the short estimation window
#' @param long_window an integer, the length of the long estimation window
#' @param panel_widths vector with 3 elements specifying the width of the
#'   vertical plot panels - parameter values (metadata), trajectories and the
#'   coverage panel
#' @return a patchwork plot
compose_coverage_patches <- function(
  plot_panels,
  short_window,
  long_window,
  panel_widths = c(0.8, 3, 2)
) {
  # We will have as many rows as the simulation scenarios
  n_rows <- length(plot_panels)

  # Set the first row as a ggplot with text only and glue the plots under each
  # other into vertical strips.
  p_trajectories <- ggplot() +
    geom_text(aes(x = 1, y = 1, label = "Simulation trajectories"), size = 5) +
    # Must specify the labels as an empty string, otherwise a later call of
    # theme() will override the labels as default "x" and "y"
    labs(x = "", y = "") +
    theme_void() +
    map(plot_panels, "trajectories") +
    plot_layout(heights = c(1, rep(6, n_rows)), axes = "collect")
  p_meta <- plot_spacer() +
    map(plot_panels, "meta") +
    plot_layout(heights = c(1, rep(6, n_rows)))

  # Compose the right panel corresponding to the empirical vs. nominal coverage
  p_coverage <- compose_subplot_by_windows(
    map(plot_panels, "coverage"),
    short_window,
    long_window
  ) +
    theme(axis.title = element_text(size = 16))

  # Collect the guides
  p_coverage <- p_coverage +
    plot_layout(
      nrow = n_rows + 1,
      ncol = 2,
      heights = c(1, rep(6, n_rows)),
      guides = "collect",
      axes = "collect"
    ) &
    get_legend_theme()

  # Create the final plot
  p_final <- (p_meta | p_trajectories | p_coverage) +
    plot_layout(guides = "collect", widths = panel_widths)
  p_final
}

#' Compose plot patches into two columns by the window length
#'
#' @description This function composes those plot patches, that have 2 versions
#'   based on the window length. Short window is on the left, long window on the
#'   right.
#' @param subplot_panels a list of lists of patchwork patches. The number of
#'   simulation scenarios is the length of the outer list. Each inner list
#'   contains 2 plots - versions of the same plot, once with the short and once
#'   with the long estimation window. The inner list elements must be
#'   accordingly named as "short" and "long".
#' @param short_window an integer, the length of the short estimation window
#' @param long_window an integer, the length of the long estimation window
#' @return a patchwork plot
compose_subplot_by_windows <- function(
  subplot_panels,
  short_window,
  long_window
) {
  n_rows <- length(subplot_panels)
  # Glue the plots under each other into vertical strips. The patchwork plot is
  # composed of two columns containing the plots of desired quantities for both
  # estimation windows in order to allow or collecting the guides and axes.
  p <- ggplot() +
    geom_text(
      aes(x = 1, y = 1, label = paste0(short_window, "-day window")),
      size = 5
    ) +
    # Must specify the labels as an empty string, otherwise a later call of
    # theme() will override the labels as default "x" and "y"
    labs(x = "", y = "") +
    theme_void() +
    map(subplot_panels, "short") +
    ggplot() +
    geom_text(
      aes(x = 1, y = 1, label = paste0(long_window, "-day window")),
      size = 5
    ) +
    # Must specify the labels as an empty string, otherwise a later call of
    # theme() will override the labels as default "x" and "y"
    labs(x = "", y = "") +
    theme_void() +
    map(subplot_panels, "long") +
    plot_layout(
      nrow = n_rows + 1,
      byrow = FALSE,
      ncol = 2,
      heights = c(1, rep(6, n_rows)),
      guides = "collect",
      axes = "collect"
    ) &
    get_legend_theme() &
    theme(axis.title = element_text(size = 16))
  p
}

#' Compose density plot patches into the final plot
#'
#' @description This function composes density plots from all the scenario runs
#'   together, including titles and labels, using the patchwork approach.
#' @param plot_panels a list of lists of patchwork patches. The number of
#'   simulation scenarios is the length of the outer list. Each inner list
#'   contains 2 lists - the density plot of R estimate and the density plot of
#'   its SE. Each list than contains one plot for the short and one plot for the
#'   long estimation window.
#' @param plot_meta_panels a list containing plots of the parameter values used
#'   in the simulation scenario.
#' @param short_window an integer, the length of the short estimation window
#' @param long_window an integer, the length of the long estimation window
#' @return a patchwork plot
compose_dens_patches <- function(
  plot_panels,
  plot_meta_panels,
  short_window,
  long_window
) {
  # Compose the density plots of R estimates
  p_R_hat <- compose_subplot_by_windows(
    plot_panels$R_hat,
    short_window,
    long_window
  )
  # Compose the density plots of SEs of R estimates
  p_se_hat <- compose_subplot_by_windows(
    plot_panels$se_hat,
    short_window,
    long_window
  )

  # Compose the panel of the metadata about the parameter values
  p_meta <- plot_spacer() + plot_meta_panels +
    plot_layout(heights = c(1, rep(6, length(plot_meta_panels))))

  # Extract the legend and remove it from the plots
  p_legend <- wrap_elements(
    ggpubr::get_legend(plot_panels$R_hat[[1]]$short & get_legend_theme())
  )
  p_R_hat <- wrap_elements(
    (p_meta | p_R_hat) +
      plot_layout(widths = c(2, 6)) &
      theme(legend.position = "none")
  )
  p_se_hat <- wrap_elements(
    (plot_spacer() | p_se_hat) +
      plot_layout(widths = c(0.25, 6)) & theme(legend.position = "none")
  )

  # Create the upper titles
  p_headline_R_hat <- plot_spacer() +
    ggplot() +
    geom_text(
      aes(x = 1, y = 1, label = "Distribution~of~hat(R)"),
      size = 6,
      parse = TRUE
    ) +
    # Must specify the labels as an empty string, otherwise a later call of
    # theme() will override the labels as default "x" and "y"
    labs(x = "", y = "") +
    theme_void() +
    plot_layout(widths = c(2, 6))
  p_headline_se_hat <- plot_spacer() +
    ggplot() +
    geom_text(
      aes(x = 1, y = 1, label = "Distribution~of~widehat(se)(hat(R))"),
      size = 6,
      parse = TRUE
    ) +
    # Must specify the labels as an empty string, otherwise a later call of
    # theme() will override the labels as default "x" and "y"
    labs(x = "", y = "") +
    theme_void() +
    plot_layout(widths = c(0.25, 6))

  # Combine into 2 big panels
  p_headline <- (p_headline_R_hat | p_headline_se_hat) +
    plot_layout(widths = c(6, 5))
  p_plot <- (p_R_hat | p_se_hat) +
    plot_layout(widths = c(6, 5))

  # Combine everything
  p_main <- wrap_elements(
    ((p_headline / p_plot) + plot_layout(heights = c(1, 22))) &
      theme(axis.title = element_text(size = 16))
  ) +
    p_legend +
    plot_layout(widths = c(10, 1))
  p_main
}

#' Compose plot patches of the overdispersion parameter estimates
#'
#' @description This function composes plots of the distribution of the
#'   overdispersion parameter estimates from all the scenario runs and also
#'   all applicable scenario blocks: NegBin-L, NegBin-L with weekday effects
#'   and NegBin-Q. The patchwork approach is used to add labels and titles.
#' @param plot_panels a list of lists of patchwork patches.
#' @param plot_meta_panels a list containing plots of the parameter values used
#'   in the simulation scenario.
#' @return a patchwork plot
compose_overdisp_patches <- function(
  plot_panels,
  plot_meta_panels
) {
  # How many vertical strips we have
  n_strips <- length(plot_panels)
  # How many rows there are in each strip
  n_rows <- length(plot_panels[[1]])

  strips <- vector("list", n_strips)
  names(strips) <- names(plot_panels)

  p_negbin_l <- plot_spacer() +
    plot_meta_panels$NegBin.L +
    # NegBin-L true data generating process
    ggplot() +
    geom_text(
      aes(x = 1, y = 1, label = "NegBin-L ground truth"),
      size = 5
    ) +
    # Must specify the labels as an empty string, otherwise a later call of
    # theme() will override the labels as default "x" and "y"
    labs(x = "", y = "") +
    theme_void() +
    plot_panels$NegBin.L_weekday_no +
    plot_layout(
      nrow = n_rows + 1,
      ncol = 2,
      byrow = FALSE,
      heights = c(1, rep(6, n_rows)),
      widths = c(2, 6),
      guides = "collect",
      axes = "collect_x",
      axis_titles = "collect"
    ) &
    theme(
      legend.position = "none"
    )

  p_negbin_q <- plot_spacer() +
    plot_meta_panels$NegBin.Q +
    # NegBin-Q true data generating process
    ggplot() +
    geom_text(
      aes(x = 1, y = 1, label = "NegBin-Q ground truth"),
      size = 5
    ) +
    labs(x = "", y = "") +
    theme_void() +
    plot_panels$NegBin.Q_weekday_no +
    plot_layout(
      nrow = n_rows + 1,
      ncol = 2,
      byrow = FALSE,
      heights = c(1, rep(6, n_rows)),
      widths = c(2, 6),
      guides = "collect",
      axes = "collect"
    ) &
    theme(
      legend.position = "none"
    )

  # Extract the legend from one of the partial plots
  p_legend <- wrap_elements(
    ggpubr::get_legend(plot_panels$NegBin.L_weekday_no[[1]])
  )

  # Combine everything
  p_main <- wrap_elements(
    ((p_negbin_l | p_negbin_q) + plot_layout(widths = c(8, 8))) &
      theme(axis.title = element_text(size = 16))
  ) +
    p_legend +
    plot_layout(widths = c(10, 1))
  p_main
}
