plot_theme <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title.position = "plot",
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(linewidth = 0.25, colour = "grey88")
    )
}

ratio_breaks <- c(0.03, 0.05, 0.1, 0.2, 0.5, 1, 2, 3, 10, 30)

forest_plot <- function(data, x_label, x_limits = NULL) {
  data <- data.table::copy(data)
  data[, label := factor(label, levels = rev(unique(label)))]
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = ratio, y = label)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dotted", linewidth = 0.45) +
    ggplot2::geom_segment(
      ggplot2::aes(x = low, xend = high, yend = label),
      linewidth = 0.55
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_log10(
      breaks = ratio_breaks,
      labels = scales::label_number(accuracy = 0.01, trim = TRUE),
      limits = x_limits
    ) +
    ggplot2::labs(x = x_label, y = NULL) +
    plot_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(hjust = 1))
  plot
}

save_plot <- function(plot, data, name, sources, kind, width, height, manifest) {
  paths <- result_paths()
  pdf <- file.path(paths$figures, paste0(name, ".pdf"))
  png <- file.path(paths$figures, paste0(name, ".png"))
  csv <- file.path(paths$figures, paste0(name, ".csv"))
  ggplot2::ggsave(pdf, plot, width = width, height = height, device = grDevices::pdf, bg = "white")
  ggplot2::ggsave(png, plot, width = width, height = height, dpi = 300, bg = "white")
  write_table(data, csv)
  source_hashes <- stats::setNames(
    vapply(sources, function(path) {
      digest::digest(file = file.path(REPOSITORY, path), algo = "sha256", serialize = FALSE)
    }, character(1)),
    sources
  )
  manifest[[name]] <- list(
    rows = nrow(data),
    data_sha256 = digest::digest(file = csv, algo = "sha256", serialize = FALSE),
    source_kind = kind,
    sources = as.list(source_hashes),
    figure_pdf = basename(pdf),
    figure_png = basename(png)
  )
  manifest
}

figure_design <- function() {
  clock <- data.table::data.table(
    regime = factor(c("Standard time", "Daylight-saving time"),
      levels = c("Daylight-saving time", "Standard time")),
    sunset = c(18, 19),
    band_start = c(18 + 10 / 60, 18 + 10 / 60),
    band_end = c(18 + 50 / 60, 18 + 50 / 60)
  )
  panel_a <- ggplot2::ggplot(clock, ggplot2::aes(y = regime)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = band_start, xend = band_end, yend = regime),
      linewidth = 7,
      alpha = 0.35
    ) +
    ggplot2::geom_point(ggplot2::aes(x = sunset), size = 2.6, shape = 18) +
    ggplot2::geom_vline(xintercept = 18, linetype = "dashed", linewidth = 0.4) +
    ggplot2::scale_x_continuous(
      breaks = c(18, 18 + 10 / 60, 18 + 50 / 60, 19),
      labels = c("18:00", "18:10", "18:50", "19:00"),
      limits = c(17.9, 19.08)
    ) +
    ggplot2::labs(
      title = "A. Illustrative solar-time exchange",
      subtitle = "Shaded civil-clock band: 18:10-18:50",
      x = "Civil-clock time",
      y = NULL
    ) +
    plot_theme() +
    ggplot2::theme(legend.position = "none")

  calendar <- data.table::data.table(
    center = factor(rep(c("Early sham", "Actual switch", "Late sham"), each = 2),
      levels = rev(c("Early sham", "Actual switch", "Late sham"))),
    period = rep(c("Before", "After"), 3),
    start = c(-56, -28, -21, 7, 14, 42),
    end = c(-43, -15, -8, 20, 27, 55)
  )
  panel_b <- ggplot2::ggplot(calendar, ggplot2::aes(y = center, colour = period)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = start, xend = end, yend = center),
      linewidth = 5,
      lineend = "butt"
    ) +
    ggplot2::geom_vline(xintercept = c(-35, 0, 35), linetype = "dotted", linewidth = 0.35) +
    ggplot2::scale_colour_manual(values = c("Before" = "black", "After" = "grey55")) +
    ggplot2::scale_x_continuous(breaks = c(-56, -35, -21, 0, 20, 35, 55)) +
    ggplot2::labs(
      title = "B. Real and symmetric sham calendars",
      x = "Days relative to the actual switch",
      y = NULL
    ) +
    plot_theme()
  list(
    plot = panel_a / panel_b + patchwork::plot_layout(heights = c(1, 1.15)),
    data = data.table::rbindlist(list(
      clock[, .(panel = "clock", label = as.character(regime), period = NA_character_, start = band_start, end = band_end, marker = sunset)],
      calendar[, .(panel = "calendar", label = as.character(center), period, start, end, marker = NA_real_)]
    ), fill = TRUE)
  )
}

figure_dynamic <- function() {
  data <- read_table("dynamic_estimates.csv", file.path(RESULTS_DIR, "estimates"))
  data[, series := paste(tools::toTitleCase(band), tools::toTitleCase(season), sep = ": ")]
  data[, week_factor := factor(week, levels = c(-4, -1, 0, 1, 2, 3))]
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = week_factor, y = ratio, group = series, shape = series)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dotted", linewidth = 0.45) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = low, ymax = high), position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.6), size = 1.8) +
    ggplot2::scale_y_log10(breaks = ratio_breaks, labels = scales::label_number(accuracy = 0.01)) +
    ggplot2::labs(
      x = "Relative week at each real or sham center",
      y = "Gain-oriented ratio of count ratios (log scale)"
    ) +
    plot_theme()
  list(plot = plot, data = data)
}

figure_primary <- function() {
  estimates <- read_table("national_models.csv", file.path(RESULTS_DIR, "estimates"))
  direct <- read_table("direct_contrasts.csv", file.path(RESULTS_DIR, "estimates"))
  rows <- list()
  for (band in BANDS) {
    rows[[length(rows) + 1L]] <- estimates[
      specification == paste0("primary_", band),
      .(label = paste(tools::toTitleCase(band), ": Pedestrian (RCR)"), ratio, low, high)
    ]
    rows[[length(rows) + 1L]] <- estimates[
      specification == paste0("other_fatal_crashes_", band),
      .(label = paste(tools::toTitleCase(band), ": Other fatal (RCR)"), ratio, low, high)
    ]
  }
  for (band in BANDS) {
    rows[[length(rows) + 1L]] <- direct[
      contrast == paste0("pedestrian_vs_other_", band),
      .(label = paste(tools::toTitleCase(band), ": Pedestrian / other"), ratio, low, high)
    ]
  }
  data <- data.table::rbindlist(rows)
  list(
    plot = forest_plot(data, "RCR, or direct ratio of RCRs (log scale)", c(0.04, 2)),
    data = data
  )
}

figure_estimator <- function() {
  primary <- read_table("national_models.csv", file.path(RESULTS_DIR, "estimates"))
  conditional <- read_table("conditional_blocks.csv", file.path(RESULTS_DIR, "estimates"))
  rows <- list()
  for (band in BANDS) {
    band_value <- band
    rows[[length(rows) + 1L]] <- primary[
      specification == paste0("primary_", band_value),
      .(label = paste(tools::toTitleCase(band_value), ": Main allocation / CR1"), ratio, low, high)
    ]
    rows[[length(rows) + 1L]] <- conditional[
      band == band_value & block_size_years == 1,
      .(label = paste(tools::toTitleCase(band_value), ": Fixed margins / year jackknife"), ratio, low, high)
    ]
    rows[[length(rows) + 1L]] <- conditional[
      band == band_value & block_size_years == 4,
      .(label = paste(tools::toTitleCase(band_value), ": Fixed margins / 4-year blocks"), ratio, low, high)
    ]
  }
  data <- data.table::rbindlist(rows)
  list(plot = forest_plot(data, "Ratio of count ratios (log scale)", c(0.05, 1.2)), data = data)
}

figure_traffic <- function() {
  traffic <- read_table("traffic_support.csv", file.path(RESULTS_DIR, "estimates"))
  labels <- data.table::data.table(
    minimum_sites = c(1, 2, 20, 30),
    label = c(
      "All matched: 44 jurisdictions / 3,892 sites",
      "≥2 or ≥10 sites: 42 jurisdictions / 3,890 sites",
      "≥20 sites: 40 jurisdictions / 3,865 sites",
      "≥30 sites: 38 jurisdictions / 3,821 sites"
    )
  )
  data <- merge(
    traffic[metric == "log_day_level" & band == "evening"],
    labels,
    by = "minimum_sites"
  )
  data <- data[, .(
    label,
    ratio = ratio_equal_jurisdiction,
    low,
    high,
    minimum_sites
  )]
  list(
    plot = forest_plot(data, "2020/2019 ratio of post/pre vehicle-volume changes", c(0.38, 1.06)),
    data = data
  )
}

figure_flow <- function() {
  flow <- support_table("source_flow.csv")
  support <- support_table("model_support.csv")[endpoint == PEDESTRIAN]
  nodes <- data.table::data.table(
    x = c(1, 1, 1, 0.55, 1.45),
    y = c(5, 4, 3, 1.7, 1.7),
    label = c(
      "National crash master\n843,530 crashes (2001–2024)",
      "Source-eligible records, including Arizona\n791,667 crashes",
      "Main study geography: 47 jurisdictions\n772,959 crashes",
      sprintf("Evening pedestrian model\n%s unique crashes\n%s appearances", format(support[band == "evening", unique_crashes_all_centers], big.mark = ","), format(support[band == "evening", event_appearances], big.mark = ",")),
      sprintf("Morning pedestrian model\n%s unique crashes\n%s appearances", format(support[band == "morning", unique_crashes_all_centers], big.mark = ","), format(support[band == "morning", event_appearances], big.mark = ","))
    )
  )
  edges <- data.table::data.table(
    x = c(1, 1, 1, 1), y = c(4.75, 3.75, 2.75, 2.75),
    xend = c(1, 1, 0.55, 1.45), yend = c(4.25, 3.25, 2.0, 2.0)
  )
  plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      arrow = ggplot2::arrow(length = grid::unit(0.12, "inches"))
    ) +
    ggplot2::geom_label(
      data = nodes,
      ggplot2::aes(x = x, y = y, label = label),
      linewidth = 0.25,
      lineheight = 0.95,
      size = 3.3
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 2), ylim = c(1.1, 5.5), clip = "off") +
    ggplot2::theme_void()
  list(plot = plot, data = nodes[, .(x, y, label)])
}

figure_counts <- function() {
  counts <- national_counts()[
    geo == "DST_main" &
      width_minutes == 40 &
      endpoint == PEDESTRIAN &
      sham == 0 &
      band %chin% BANDS
  ]
  data <- counts[, .(events = sum(events)), by = .(season, band, period)]
  data[, comparison := paste(tools::toTitleCase(season), tools::toTitleCase(band))]
  data[, period := factor(period, levels = c("pre", "post"), labels = c("Before", "After"))]
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = comparison, y = events, fill = period)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72), width = 0.64) +
    ggplot2::scale_fill_grey(start = 0.25, end = 0.75) +
    ggplot2::labs(x = NULL, y = "Pedestrian-fatal crash count", fill = NULL) +
    plot_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 18, hjust = 1))
  list(plot = plot, data = data)
}

figure_legacy <- function() {
  data <- read_table("legacy_weekly_estimates.csv", file.path(RESULTS_DIR, "estimates"))
  data[, series := paste(tools::toTitleCase(band), tools::toTitleCase(season), sep = ": ")]
  data[, week_factor := factor(week, levels = c(-4, -3, -2, 0, 1, 2, 3))]
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = week_factor, y = ratio, group = series, shape = series)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dotted", linewidth = 0.45) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = low, ymax = high), position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.6), size = 1.8) +
    ggplot2::scale_y_log10(breaks = ratio_breaks, labels = scales::label_number(accuracy = 0.01)) +
    ggplot2::labs(
      x = "Relative week at the actual clock change (Monday–Saturday)",
      y = "Focal/reference ratio relative to week −1 (log scale)"
    ) +
    plot_theme()
  list(plot = plot, data = data)
}

figure_leave_year <- function() {
  data <- read_table("leave_one_year.csv", file.path(RESULTS_DIR, "estimates"))
  data[, band_label := factor(
    tools::toTitleCase(band),
    levels = c("Evening", "Morning")
  )]
  bands <- unique(data[, .(band_label, full_ratio)])
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = omitted_policy_year, y = ratio)) +
    ggplot2::geom_hline(
      data = bands,
      ggplot2::aes(yintercept = full_ratio),
      linetype = "dashed",
      linewidth = 0.45
    ) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~band_label, ncol = 1, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = seq(2001, 2024, by = 2)) +
    ggplot2::labs(
      x = "Omitted policy-transition year",
      y = "Leave-one-year RCR"
    ) +
    plot_theme() +
    ggplot2::theme(legend.position = "none")
  list(plot = plot, data = data[, !c("low", "high")])
}

figure_calendar <- function() {
  data <- read_table("calendar_sensitivity.csv", file.path(RESULTS_DIR, "estimates"))
  labels <- c(
    matched_14d_gap35 = "14 days; centers ±35 days",
    nonoverlap_14d_gap42 = "14 days; centers ±42 days (no overlap)",
    short_7d_gap35 = "7 days; centers ±35 days",
    long_21d_gap35 = "21 days; centers ±35 days"
  )
  data[, label := labels[specification]]
  data[, label := factor(label, levels = rev(unname(labels)))]
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = ratio, y = label, shape = band)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dotted", linewidth = 0.45) +
    ggplot2::geom_segment(
      ggplot2::aes(x = low, xend = high, yend = label),
      position = ggplot2::position_dodge(width = 0.45),
      linewidth = 0.55
    ) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.45), size = 2) +
    ggplot2::scale_x_log10(breaks = c(0.03, 0.1, 0.3, 1), limits = c(0.03, 1.5)) +
    ggplot2::labs(x = "Gain-oriented ratio of count ratios (log scale)", y = NULL, shape = NULL) +
    plot_theme()
  list(plot = plot, data = data)
}

render_all_figures <- function() {
  manifest <- list()
  figures <- list(
    Figure_1_design = c(figure_design(), list(
      sources = c("data/calendar.csv"), kind = "declared clock example and calendar", width = 9, height = 7
    )),
    Figure_2_aligned_dynamics = c(figure_dynamic(), list(
      sources = c("results/R/estimates/dynamic_estimates.csv"), kind = "re-estimated dynamic contrasts", width = 9, height = 5.5
    )),
    Figure_3_primary_contrasts = c(figure_primary(), list(
      sources = c("results/R/estimates/national_models.csv", "results/R/estimates/direct_contrasts.csv"), kind = "re-estimated main contrasts", width = 9.2, height = 5.2
    )),
    Figure_S1_estimator_sensitivity = c(figure_estimator(), list(
      sources = c("results/R/estimates/national_models.csv", "results/R/estimates/conditional_blocks.csv"), kind = "estimator sensitivity", width = 9.2, height = 5.2
    )),
    Figure_S2_traffic_support = c(figure_traffic(), list(
      sources = c("results/R/estimates/traffic_support.csv"), kind = "traffic support sensitivity", width = 9.2, height = 4.2
    )),
    Figure_S3_sample_flow = c(figure_flow(), list(
      sources = c("data/support/source_flow.csv", "data/support/model_support.csv"), kind = "source and analytic support", width = 8.5, height = 7
    )),
    Figure_S4_seasonal_counts = c(figure_counts(), list(
      sources = c("data/national.csv.gz"), kind = "raw aggregate counts", width = 8.5, height = 4.8
    )),
    Figure_S5_original_event_weeks = c(figure_legacy(), list(
      sources = c("data/sensitivity/legacy_weeks.csv", "results/R/estimates/legacy_weekly_estimates.csv"), kind = "re-estimated legacy profiles", width = 9, height = 5.5
    )),
    Figure_S6_calendar = c(figure_calendar(), list(
      sources = c("results/R/estimates/calendar_sensitivity.csv"), kind = "calendar sensitivity", width = 9, height = 4.6
    )),
    Figure_S7_policy_year_influence = c(figure_leave_year(), list(
      sources = c("results/R/estimates/leave_one_year.csv"), kind = "leave-one-policy-year main-model estimates; dashed lines are full-model estimates", width = 9.2, height = 6.6
    ))
  )

  for (name in names(figures)) {
    item <- figures[[name]]
    manifest <- save_plot(
      item$plot,
      item$data,
      name,
      item$sources,
      item$kind,
      item$width,
      item$height,
      manifest
    )
  }
  jsonlite::write_json(
    manifest,
    file.path(RESULTS_DIR, "figures", "FIGURE_DATA_MANIFEST.json"),
    pretty = TRUE,
    auto_unbox = TRUE
  )
  invisible(manifest)
}
