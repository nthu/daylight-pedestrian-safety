source(file.path("R", "00_setup.R"))

reference_sensitivity <- function() {
  counts <- reference_counts()
  specifications <- list(
    list(name = "evening_ref_early_evening", band = "evening", reference = "early_evening"),
    list(name = "evening_ref_late_evening", band = "evening", reference = "late_evening"),
    list(name = "evening_DST_holidaypaired", band = "evening", geography = "DST_holidaypaired"),
    list(name = "evening_DST_weekday", band = "evening", geography = "DST_weekday"),
    list(name = "evening_DST_weekend", band = "evening", geography = "DST_weekend"),
    list(name = "morning_ref_early_morning", band = "morning", reference = "early_morning"),
    list(name = "morning_ref_late_morning", band = "morning", reference = "late_morning"),
    list(name = "morning_DST_holidaypaired", band = "morning", geography = "DST_holidaypaired"),
    list(name = "morning_DST_weekday", band = "morning", geography = "DST_weekday"),
    list(name = "control_night_vs_noon", band = "night", reference = "noon"),
    list(name = "control_late_vs_early_evening", band = "late_evening", reference = "early_evening")
  )
  data.table::rbindlist(lapply(specifications, function(specification) {
    name <- specification$name
    fit <- do.call(
      estimate_allocation,
      c(list(counts = counts), specification[names(specification) != "name"])
    )
    data.table::data.table(
      specification = name,
      ratio = fit$ratio,
      low = fit$low,
      high = fit$high,
      p_value = fit$p_value,
      se = fit$se,
      df = fit$df
    )
  }))
}

calendar_sensitivity <- function() {
  rows <- list()
  for (name in names(calendar_specifications)) {
    counts <- sensitivity_table(calendar_specifications[[name]])
    geography <- counts$geo[1]
    for (band in BANDS) {
      fit <- estimate_allocation(counts, band = band, geography = geography)
      rows[[length(rows) + 1L]] <- data.table::data.table(
        specification = name,
        band = band,
        ratio = fit$ratio,
        low = fit$low,
        high = fit$high,
        se = fit$se,
        df = fit$df
      )
    }
  }
  data.table::rbindlist(rows)
}

dynamic_sensitivity <- function() {
  rows <- list()
  for (week in c(-4, -1, 0, 1, 2, 3)) {
    filename <- sprintf("dynamic_counts_week_%+d.csv", week)
    counts <- sensitivity_table(filename)
    geography <- counts$geo[1]
    for (band in BANDS) {
      for (season in SEASONS) {
        fit <- estimate_allocation(
          counts,
          band = band,
          seasons = season,
          geography = geography
        )
        rows[[length(rows) + 1L]] <- data.table::data.table(
          week = week,
          band = band,
          season = season,
          ratio = fit$ratio,
          low = fit$low,
          high = fit$high,
          se = fit$se,
          df = fit$df
        )
      }
    }
  }
  data.table::rbindlist(rows)
}

whole_day_analysis <- function(counts) {
  rows <- list()
  for (endpoint_value in c(PEDESTRIAN, OTHER_FATAL, "all_fatal_crashes")) {
    selected <- counts[
      geo == "DST_main" &
        width_minutes == 40 &
        band == "all_day" &
        endpoint == endpoint_value
    ]
    table <- data.table::dcast(
      selected,
      year + season ~ sham + period,
      value.var = "events"
    )
    for (season_value in SEASONS) {
      annual <- table[season == season_value]
      required <- as.vector(outer(c(-35, 0, 35), c("pre", "post"), paste, sep = "_"))
      if (any(annual[, ..required] <= 0)) {
        stop("Whole-day counts contain zero cells.", call. = FALSE)
      }
      changes <- lapply(c(-35, 0, 35), function(center) {
        log(annual[[paste0(center, "_post")]] / annual[[paste0(center, "_pre")]])
      })
      names(changes) <- c("-35", "0", "35")
      methods <- list(
        raw_actual = changes[["0"]],
        symmetric_sham_adjusted = changes[["0"]] -
          0.5 * (changes[["-35"]] + changes[["35"]])
      )
      for (method in names(methods)) {
        values <- methods[[method]]
        standard_error <- stats::sd(values) / sqrt(length(values))
        interval <- t_interval(mean(values), standard_error, length(values) - 1L)
        rows[[length(rows) + 1L]] <- data.table::data.table(
          endpoint = endpoint_value,
          season = season_value,
          method = method,
          ratio = interval$ratio,
          low = interval$low,
          high = interval$high
        )
      }
    }
  }
  data.table::rbindlist(rows)
}

pre2022_analysis <- function(counts) {
  rows <- list()
  for (band in BANDS) {
    fitted <- list()
    for (endpoint in ENDPOINTS) {
      fit <- estimate_allocation(
        counts,
        band = band,
        endpoint = endpoint,
        years = 2001:2020
      )
      fitted[[endpoint]] <- fit
      rows[[length(rows) + 1L]] <- data.table::data.table(
        band = band,
        endpoint = endpoint,
        ratio = fit$ratio,
        low = fit$low,
        high = fit$high,
        p_value = fit$p_value,
        se = fit$se,
        df = fit$df
      )
    }
    direct <- contrast_allocation(
      fitted[[PEDESTRIAN]],
      fitted[[OTHER_FATAL]],
      paste0("pre2022_", band)
    )
    rows[[length(rows) + 1L]] <- data.table::data.table(
      band = band,
      endpoint = "direct",
      ratio = direct$ratio,
      low = direct$low,
      high = direct$high,
      p_value = direct$p_value,
      se = direct$se,
      df = direct$df
    )
  }
  data.table::rbindlist(rows)
}

clock_shift_analysis <- function() {
  counts <- sensitivity_table("clock_shift.csv.gz")
  rows <- list()
  for (minutes in sort(unique(counts$clock_shift_minutes))) {
    shifted <- counts[clock_shift_minutes == minutes]
    for (band in BANDS) {
      fitted <- list()
      for (endpoint in ENDPOINTS) {
        fit <- estimate_allocation(shifted, band = band, endpoint = endpoint)
        fitted[[endpoint]] <- fit
        rows[[length(rows) + 1L]] <- data.table::data.table(
          clock_shift_minutes = minutes,
          band = band,
          endpoint = endpoint,
          ratio = fit$ratio,
          low = fit$low,
          high = fit$high,
          se = fit$se,
          df = fit$df
        )
      }
      direct <- contrast_allocation(
        fitted[[PEDESTRIAN]],
        fitted[[OTHER_FATAL]],
        "clock_shift"
      )
      rows[[length(rows) + 1L]] <- data.table::data.table(
        clock_shift_minutes = minutes,
        band = band,
        endpoint = "direct",
        ratio = direct$ratio,
        low = direct$low,
        high = direct$high,
        se = direct$se,
        df = direct$df
      )
    }
  }
  data.table::rbindlist(rows)
}

run_sensitivity <- function(full = TRUE) {
  paths <- result_paths()
  counts <- national_counts()
  exact_blocks <- data.table::rbindlist(lapply(BANDS, function(band) {
    exact_block_intervals(
      primary_tables(),
      band,
      block_sizes = if (full) 1:4 else c(1, 4)
    )
  }))
  triplets <- triplet_analysis(counts)
  common_time <- common_time_analysis(counts)

  curvature <- five_center_analysis(full = full)
  diagnostics <- review_diagnostics(counts)
  outputs <- list(
    reference_sensitivity = reference_sensitivity(),
    calendar_sensitivity = calendar_sensitivity(),
    dynamic_estimates = dynamic_sensitivity(),
    conditional_blocks = exact_blocks,
    whole_day = whole_day_analysis(counts),
    clock_shift = clock_shift_analysis(),
    wild_cluster_calibration = wild_cluster_calibration(counts),
    five_center_curvature = curvature$five_center_curvature,
    five_center_deletions = curvature$five_center_deletions,
    season_heterogeneity = diagnostics$season_heterogeneity,
    leave_one_year = diagnostics$leave_one_year,
    policy_era_interaction = diagnostics$policy_era_interaction,
    period_interactions = diagnostics$period_interactions,
    arizona_negative_control = diagnostics$arizona_negative_control,
    alaska_inclusion = diagnostics$alaska_inclusion,
    pre2022 = pre2022_analysis(counts),
    annual_slope_conditional = triplets$estimates,
    triplet_support = triplets$support,
    common_time = common_time$estimates,
    sham_decomposition = common_time$decomposition,
    standardization_support = common_time$support
  )
  for (name in names(outputs)) {
    write_table(outputs[[name]], file.path(paths$estimates, paste0(name, ".csv")))
  }

  cat("\nKey calendar sensitivities\n")
  print(outputs$calendar_sensitivity)
  cat("\nCommon-time direct comparisons\n")
  print(outputs$common_time[endpoint == "direct", .(band, ratio, low, high)])
  invisible(outputs)
}

if (sys.nframe() == 0L) run_sensitivity(full = TRUE)
