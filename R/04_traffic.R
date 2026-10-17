source(file.path("R", "00_setup.R"))

run_traffic <- function() {
  traffic <- traffic_counts()
  rows <- list()
  metrics <- c("log_day_level", "log_focal_daily_share", "log_focal_vs_midnight")
  for (metric_value in metrics) {
    for (band_value in BANDS) {
      for (minimum_sites in c(1, 2, 5, 10, 20, 30)) {
        selected <- traffic[
          metric == metric_value &
            band == band_value &
            n_sites >= minimum_sites
        ]
        jurisdictions <- nrow(selected)
        if (jurisdictions < 2L) next
        values <- selected$shock_log_ratio_2020_vs_2019
        mean_value <- mean(values)
        standard_error <- stats::sd(values) / sqrt(jurisdictions)
        interval <- t_interval(mean_value, standard_error, jurisdictions - 1L)
        rows[[length(rows) + 1L]] <- data.table::data.table(
          metric = metric_value,
          band = band_value,
          minimum_sites = minimum_sites,
          jurisdictions = jurisdictions,
          sites = sum(selected$n_sites),
          ratio_equal_jurisdiction = interval$ratio,
          low = interval$low,
          high = interval$high,
          ratio_median_jurisdiction = exp(stats::median(values)),
          ratio_site_weighted = exp(stats::weighted.mean(values, selected$n_sites))
        )
      }
    }
  }
  results <- data.table::rbindlist(rows)
  paths <- result_paths()
  write_table(results, file.path(paths$estimates, "traffic_support.csv"))

  cat("\nMatched-traffic sensitivity: daily volume, evening\n")
  print(results[
    metric == "log_day_level" & band == "evening",
    .(minimum_sites, jurisdictions, sites, ratio_equal_jurisdiction, low, high)
  ])
  invisible(results)
}

if (sys.nframe() == 0L) run_traffic()
