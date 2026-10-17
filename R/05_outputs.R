source(file.path("R", "00_setup.R"))

run_outputs <- function() {
  paths <- result_paths()
  legacy <- estimate_legacy_weeks()
  write_table(
    legacy$estimates,
    file.path(paths$estimates, "legacy_weekly_estimates.csv")
  )
  write_table(
    legacy$pretests,
    file.path(paths$estimates, "legacy_weekly_pretests.csv")
  )

  national <- read_table("national_models.csv", paths$estimates)
  direct <- read_table("direct_contrasts.csv", paths$estimates)
  table_rows <- list()
  endpoints <- list(
    "Fatal pedestrian crashes" = "primary_",
    "Other fatal crashes" = "other_fatal_crashes_",
    "All fatal crashes" = "all_fatal_crashes_"
  )
  for (label in names(endpoints)) {
    prefix <- endpoints[[label]]
    table_rows[[length(table_rows) + 1L]] <- data.table::data.table(
      endpoint = label,
      evening = national[specification == paste0(prefix, "evening"),
        sprintf("%.3f (%.3f–%.3f)", ratio, low, high)],
      morning = national[specification == paste0(prefix, "morning"),
        sprintf("%.3f (%.3f–%.3f)", ratio, low, high)]
    )
  }
  table_rows[[length(table_rows) + 1L]] <- data.table::data.table(
    endpoint = "Pedestrian versus other: ratio of RCRs",
    evening = direct[contrast == "pedestrian_vs_other_evening",
      sprintf("%.3f (%.3f–%.3f)", ratio, low, high)],
    morning = direct[contrast == "pedestrian_vs_other_morning",
      sprintf("%.3f (%.3f–%.3f)", ratio, low, high)]
  )
  table2 <- data.table::rbindlist(table_rows)
  write_table(table2, file.path(paths$tables, "table_2_primary_contrasts.csv"))

  leave_year <- read_table("leave_one_year.csv", paths$estimates)
  conditional <- read_table("conditional_blocks.csv", paths$estimates)
  wild <- read_table("wild_cluster_calibration.csv", paths$estimates)
  five <- read_table("five_center_curvature.csv", paths$estimates)
  calendar <- read_table("calendar_sensitivity.csv", paths$estimates)

  jackknife_values <- function(band_value) {
    values <- log(leave_year[band == band_value, ratio])
    blocks <- length(values)
    standard_error <- sqrt((blocks - 1) / blocks * sum((values - mean(values))^2))
    beta <- log(national[specification == paste0("primary_", band_value), ratio])
    critical <- stats::qt(0.975, blocks - 1)
    c(exp(beta), exp(beta - critical * standard_error), exp(beta + critical * standard_error))
  }

  specifications <- list(
    list(label = "Main three-center allocation, CR1", source = "main", scope = "2001–2024"),
    list(label = "Exact policy-year multiplier calibration", source = "wild", scope = "2001–2024"),
    list(label = "Main allocation, policy-year jackknife", source = "jackknife", scope = "2001–2024"),
    list(label = "Fixed-margin conditional, policy-year jackknife", source = "conditional", scope = "2001–2024"),
    list(label = "Fixed comparison-center solar threshold", source = "fixed", scope = "2001–2024"),
    list(label = "Nonoverlapping ±42-day calendar", source = "nonoverlap", scope = "2002–2023"),
    list(label = "Five-center quadratic background", source = "five", scope = "2002–2023")
  )
  diagnostic_rows <- lapply(specifications, function(specification) {
    row <- list(analysis = specification$label, scope = specification$scope)
    for (band_value in BANDS) {
      values <- switch(
        specification$source,
        main = national[specification == paste0("primary_", band_value), c(ratio, low, high)],
        wild = wild[comparison == paste0("pedestrian_", band_value), c(ratio, wild_low, wild_high)],
        jackknife = jackknife_values(band_value),
        conditional = conditional[band == band_value & block_size_years == 1, c(ratio, low, high)],
        fixed = national[specification == paste0("fixed_clock_", band_value), c(ratio, low, high)],
        nonoverlap = calendar[specification == "nonoverlap_14d_gap42" & band == band_value, c(ratio, low, high)],
        five = five[model == "five_center_quadratic" & endpoint == PEDESTRIAN & band == band_value, c(ratio, low, high)]
      )
      values <- as.numeric(values)
      row[[paste0(band_value, "_ratio")]] <- values[1]
      row[[paste0(band_value, "_low")]] <- values[2]
      row[[paste0(band_value, "_high")]] <- values[3]
    }
    data.table::as.data.table(row)
  })
  table3 <- data.table::rbindlist(diagnostic_rows)
  write_table(table3, file.path(paths$tables, "table_3_key_diagnostics.csv"))

  manifest <- render_all_figures()
  cat("\nGenerated", length(manifest), "figures in", paths$figures, "\n")
  invisible(list(table2 = table2, figures = manifest))
}

if (sys.nframe() == 0L) run_outputs()
