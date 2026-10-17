mantel_haenszel_log_contrast <- function(
    frame,
    unit,
    band,
    endpoint,
    years) {
  band_value <- band
  endpoint_value <- endpoint
  selected <- frame[
    band %chin% c(band_value, "night") &
      endpoint == endpoint_value &
      year %in% years
  ]
  formula <- stats::as.formula(
    paste(unit, "+ year + season + sham ~ period + band")
  )
  table <- data.table::dcast(selected, formula, value.var = "events")
  required <- c(
    paste0("post_", band),
    "post_night",
    paste0("pre_", band),
    "pre_night"
  )
  if (!all(required %in% names(table))) {
    stop("Fixed-unit table is incomplete.", call. = FALSE)
  }

  post_focal <- table[[paste0("post_", band)]]
  post_reference <- table$post_night
  pre_focal <- table[[paste0("pre_", band)]]
  pre_reference <- table$pre_night
  total <- post_focal + post_reference + pre_focal + pre_reference
  numerator <- fifelse(total > 0, post_focal * pre_reference / total, 0)
  denominator <- fifelse(total > 0, post_reference * pre_focal / total, 0)
  table[, `:=`(numerator = numerator, denominator = denominator)]
  sums <- table[, .(
    numerator = sum(numerator),
    denominator = sum(denominator)
  ), by = .(season, sham)]
  if (any(sums$numerator <= 0) || any(sums$denominator <= 0)) {
    stop("Fixed-unit aggregate support is zero.", call. = FALSE)
  }
  sums[, log_odds := log(numerator / denominator)]

  seasonal <- numeric()
  for (season_value in SEASONS) {
    direction <- ifelse(season_value == "spring", 1, -1) *
      ifelse(band == "evening", 1, -1)
    values <- sums[season == season_value]
    data.table::setkey(values, sham)
    seasonal <- c(
      seasonal,
      direction * (
        values[J(0), log_odds] -
          0.5 * (values[J(-35), log_odds] + values[J(35), log_odds])
      )
    )
  }
  list(log_ratio = mean(seasonal), table = table)
}

information_support <- function(table, band, unit) {
  a <- table[[paste0("post_", band)]]
  b <- table$post_night
  c0 <- table[[paste0("pre_", band)]]
  d <- table$pre_night
  total <- a + b + c0 + d
  post_total <- a + b
  focal_total <- a + c0
  informative <- pmax(0, post_total - (total - focal_total)) <
    pmin(post_total, focal_total)
  table[, c(unit, "year", "season", "sham"), with = FALSE][, `:=`(
    appearances = total,
    informative = informative,
    all_zero = total == 0
  )]
}

state_analysis <- function(counts) {
  years <- sort(unique(counts$year))
  endpoints <- c(PEDESTRIAN, OTHER_FATAL, "all_fatal_crashes")
  estimates <- list()
  deletions <- list()
  supports <- list()
  stored <- list()

  for (endpoint in endpoints) {
    for (band in BANDS) {
      full <- mantel_haenszel_log_contrast(
        counts,
        unit = "state",
        band = band,
        endpoint = endpoint,
        years = years
      )
      deleted <- vapply(
        years,
        function(year) mantel_haenszel_log_contrast(
          counts,
          unit = "state",
          band = band,
          endpoint = endpoint,
          years = setdiff(years, year)
        )$log_ratio,
        numeric(1)
      )
      interval <- jackknife_interval(full$log_ratio, deleted)
      estimates[[length(estimates) + 1L]] <- data.table::data.table(
        endpoint = endpoint,
        band = band,
        ratio = interval$ratio,
        log_ratio = interval$log_ratio,
        se = interval$se,
        low = interval$low,
        high = interval$high,
        p_value = interval$p_value,
        df = interval$df,
        states = 47L
      )
      stored[[paste(endpoint, band, sep = "|")]] <- list(
        log_ratio = full$log_ratio,
        deleted = deleted
      )
      deletions[[length(deletions) + 1L]] <- data.table::data.table(
        endpoint = endpoint,
        band = band,
        omitted_policy_year = years,
        log_ratio = deleted
      )
      support <- information_support(full$table, band, "state")
      support[, `:=`(endpoint = endpoint, band = band)]
      supports[[length(supports) + 1L]] <- support
    }
  }

  direct <- list()
  for (band in BANDS) {
    pedestrian <- stored[[paste(PEDESTRIAN, band, sep = "|")]]
    other <- stored[[paste(OTHER_FATAL, band, sep = "|")]]
    interval <- jackknife_interval(
      pedestrian$log_ratio - other$log_ratio,
      pedestrian$deleted - other$deleted
    )
    direct[[length(direct) + 1L]] <- data.table::data.table(
      band = band,
      ratio = interval$ratio,
      log_ratio = interval$log_ratio,
      se = interval$se,
      low = interval$low,
      high = interval$high,
      p_value = interval$p_value,
      df = interval$df
    )
  }

  support <- data.table::rbindlist(supports)
  summary <- support[, .(
    appearances = sum(appearances),
    tables = .N,
    informative_tables = sum(informative),
    all_zero_tables = sum(all_zero)
  ), by = .(endpoint, band, state)]
  list(
    state_effects = data.table::rbindlist(estimates),
    state_direct = data.table::rbindlist(direct),
    state_deletions = data.table::rbindlist(deletions),
    state_support = summary
  )
}

mississippi_analysis <- function(counts) {
  years <- sort(unique(counts$year))
  endpoints <- c(
    "all_crashes",
    "injury_or_fatal_crashes",
    "nonfatal_injury_crashes",
    "property_damage_only_crashes"
  )
  estimates <- list()
  for (endpoint in endpoints) {
    for (band in BANDS) {
      full <- mantel_haenszel_log_contrast(
        counts,
        unit = "county_key",
        band = band,
        endpoint = endpoint,
        years = years
      )
      deleted <- vapply(
        years,
        function(year) mantel_haenszel_log_contrast(
          counts,
          unit = "county_key",
          band = band,
          endpoint = endpoint,
          years = setdiff(years, year)
        )$log_ratio,
        numeric(1)
      )
      interval <- jackknife_interval(full$log_ratio, deleted)
      estimates[[length(estimates) + 1L]] <- data.table::data.table(
        endpoint = endpoint,
        band = band,
        ratio = interval$ratio,
        log_ratio = interval$log_ratio,
        se = interval$se,
        low = interval$low,
        high = interval$high,
        p_value = interval$p_value,
        df = interval$df
      )
    }
  }
  effects <- data.table::rbindlist(estimates)
  effects[, holm_p := holm_adjust(p_value)]

  pooled <- counts[, .(events = sum(events)), by = .(
    year,
    season,
    sham,
    period,
    band,
    width_minutes,
    endpoint
  )]
  pooled[, geo := "MS_fixed_county"]

  direct <- list()
  for (band in BANDS) {
    injury <- estimate_allocation(
      pooled,
      band = band,
      endpoint = "nonfatal_injury_crashes",
      geography = "MS_fixed_county"
    )
    property_only <- estimate_allocation(
      pooled,
      band = band,
      endpoint = "property_damage_only_crashes",
      geography = "MS_fixed_county"
    )
    comparison <- contrast_allocation(
      injury,
      property_only,
      "nonfatal_injury_vs_property_damage_only"
    )
    direct[[length(direct) + 1L]] <- data.table::as.data.table(
      c(list(band = band), comparison)
    )
  }

  variance <- list()
  for (band in BANDS) {
    full <- estimate_allocation(
      pooled,
      band = band,
      endpoint = "all_crashes",
      geography = "MS_fixed_county"
    )
    deleted <- vapply(
      years,
      function(year) estimate_allocation(
        pooled,
        band = band,
        endpoint = "all_crashes",
        geography = "MS_fixed_county",
        years = setdiff(years, year)
      )$log_ratio,
      numeric(1)
    )
    jackknife <- jackknife_interval(full$log_ratio, deleted)
    variance[[length(variance) + 1L]] <- data.table::data.table(
      band = band,
      method = "allocation_CR1",
      ratio = full$ratio,
      se = full$se,
      low = full$low,
      high = full$high
    )
    variance[[length(variance) + 1L]] <- data.table::data.table(
      band = band,
      method = "allocation_year_jackknife",
      ratio = jackknife$ratio,
      se = jackknife$se,
      low = jackknife$low,
      high = jackknife$high
    )
  }

  list(
    mississippi_effects = effects,
    mississippi_direct = data.table::rbindlist(direct, fill = TRUE),
    mississippi_variance = data.table::rbindlist(variance)
  )
}
