seasonal_heterogeneity <- function(counts = national_counts()) {
  rows <- list()
  for (band_value in BANDS) {
    fitted <- list()
    for (season_value in c("spring", "autumn")) {
      result <- estimate_allocation(
        counts,
        band = band_value,
        seasons = season_value
      )
      fitted[[season_value]] <- result
      rows[[length(rows) + 1L]] <- data.table::data.table(
        band = band_value,
        result_type = "season_estimate",
        label = season_value,
        ratio = result$ratio,
        low = result$low,
        high = result$high,
        se = result$se,
        p_value = result$p_value,
        df = result$df,
        shared_years = length(result$years)
      )
    }
    comparison <- contrast_allocation(
      fitted$spring,
      fitted$autumn,
      paste0("spring_vs_autumn_", band_value)
    )
    rows[[length(rows) + 1L]] <- data.table::data.table(
      band = band_value,
      result_type = "season_contrast",
      label = "spring_vs_autumn",
      ratio = comparison$ratio,
      low = comparison$low,
      high = comparison$high,
      se = comparison$se,
      p_value = comparison$p_value,
      df = comparison$df,
      shared_years = comparison$shared_years
    )
  }
  data.table::rbindlist(rows)
}

leave_one_year_influence <- function(counts = national_counts()) {
  rows <- list()
  years <- 2001:2024
  for (band_value in BANDS) {
    full <- estimate_allocation(counts, band = band_value)
    for (omitted in years) {
      fitted <- estimate_allocation(
        counts,
        band = band_value,
        years = setdiff(years, omitted)
      )
      rows[[length(rows) + 1L]] <- data.table::data.table(
        band = band_value,
        omitted_policy_year = omitted,
        ratio = fitted$ratio,
        low = fitted$low,
        high = fitted$high,
        full_ratio = full$ratio,
        full_low = full$low,
        full_high = full$high
      )
    }
  }
  data.table::rbindlist(rows)
}

fit_era_model <- function(counts, band, groups) {
  band_value <- band
  selected <- counts[
    geo == "DST_main" &
      width_minutes == 40 &
      endpoint == PEDESTRIAN &
      band %chin% c(band_value, "night")
  ]
  table <- data.table::dcast(
    selected,
    year + season + sham ~ period + band,
    value.var = "events"
  )
  data.table::setorderv(table, c("year", "season", "sham"))

  a <- table[[paste0("post_", band)]]
  b <- table$post_night
  c0 <- table[[paste0("pre_", band)]]
  d <- table$pre_night
  total_2x2 <- a + b + c0 + d
  focal_total <- a + c0
  post_total <- a + b
  informative <- pmax(0, post_total - (total_2x2 - focal_total)) < pmin(post_total, focal_total)
  table <- table[informative]
  strata <- nrow(table)

  columns <- list()
  extra_names <- character()
  present <- unique(table[, .(year, season)])
  data.table::setorderv(present, c("year", "season"))
  for (index in seq_len(nrow(present))) {
    year_value <- present$year[index]
    season_value <- present$season[index]
    columns[[length(columns) + 1L]] <- as.numeric(
      table$year == year_value & table$season == season_value
    )
    extra_names <- c(extra_names, paste0("post_", year_value, "_", season_value))
  }
  for (season_value in sort(unique(table$season))) {
    columns[[length(columns) + 1L]] <-
      as.numeric(table$season == season_value) * table$sham / 35
    extra_names <- c(extra_names, paste0("post_", season_value, "_centre_slope"))
  }
  direction <- ifelse(table$season == "spring", 1, -1) * ifelse(band == "evening", 1, -1)
  for (group in groups) {
    label <- group$label
    selector <- group$selector
    active <- as.numeric(vapply(table$year, selector, logical(1)))
    columns[[length(columns) + 1L]] <- as.numeric(table$sham == 0) * direction * active
    extra_names <- c(extra_names, paste0("effect_", label))
  }
  extras <- do.call(cbind, columns)

  design <- matrix(0, nrow = 2L * strata, ncol = strata + ncol(extras))
  success <- failure <- numeric(2L * strata)
  policy_year <- integer(2L * strata)
  row <- 1L
  for (index in seq_len(strata)) {
    for (post in 0:1) {
      design[row, index] <- 1
      if (post == 1L) design[row, (strata + 1L):ncol(design)] <- extras[index, ]
      period <- if (post == 0L) "pre" else "post"
      success[row] <- table[[paste0(period, "_", band)]][index]
      failure[row] <- table[[paste0(period, "_night")]][index]
      policy_year[row] <- table$year[index]
      row <- row + 1L
    }
  }

  names <- c(paste0("stratum_", seq_len(strata) - 1L), extra_names)
  colnames(design) <- names
  fit <- stats::glm.fit(
    x = design,
    y = cbind(success, failure),
    family = stats::binomial(),
    control = stats::glm.control(epsilon = 1e-11, maxit = 500)
  )
  if (!isTRUE(fit$converged) || any(!is.finite(fit$coefficients))) {
    stop("Era-interacted grouped-binomial model failed.", call. = FALSE)
  }

  totals <- success + failure
  probability <- fit$fitted.values
  score_rows <- design * (success - totals * probability)
  information <- crossprod(
    design,
    design * (totals * probability * (1 - probability))
  )
  eigenvalues <- eigen((information + t(information)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) <= max(eigenvalues) * 1e-12) {
    stop("Era-interacted information matrix is not full rank.", call. = FALSE)
  }
  bread <- solve(information)
  cluster_scores <- rowsum(score_rows, policy_year, reorder = TRUE)
  clusters <- as.integer(rownames(cluster_scores))
  G <- nrow(cluster_scores)
  N <- nrow(design)
  K <- ncol(design)
  correction <- G / (G - 1) * (N - 1) / (N - K)
  covariance <- bread %*% crossprod(cluster_scores) %*% bread * correction
  df <- G - 1L
  critical <- stats::qt(0.975, df)

  effects <- list()
  for (group in groups) {
    label <- group$label
    selector <- group$selector
    index <- match(paste0("effect_", label), names)
    beta <- fit$coefficients[index]
    se <- sqrt(covariance[index, index])
    effects[[length(effects) + 1L]] <- data.table::data.table(
      band = band,
      result_type = "era_coefficient",
      label = label,
      log_ratio = beta,
      ratio = exp(beta),
      se = se,
      low = exp(beta - critical * se),
      high = exp(beta + critical * se),
      p_value = 2 * stats::pt(-abs(beta / se), df),
      df = df,
      years_in_era = sum(vapply(clusters, selector, logical(1))),
      from = "",
      to = ""
    )
  }

  comparisons <- list()
  for (first in seq_len(length(groups) - 1L)) {
    for (second in (first + 1L):length(groups)) {
      label_a <- groups[[first]]$label
      label_b <- groups[[second]]$label
      index_a <- match(paste0("effect_", label_a), names)
      index_b <- match(paste0("effect_", label_b), names)
      beta <- fit$coefficients[index_b] - fit$coefficients[index_a]
      se <- sqrt(
        covariance[index_b, index_b] +
          covariance[index_a, index_a] -
          2 * covariance[index_b, index_a]
      )
      comparisons[[length(comparisons) + 1L]] <- data.table::data.table(
        band = band,
        result_type = "interaction_difference",
        label = paste0(label_b, "_vs_", label_a),
        log_ratio = beta,
        ratio = exp(beta),
        se = se,
        low = exp(beta - critical * se),
        high = exp(beta + critical * se),
        p_value = 2 * stats::pt(-abs(beta / se), df),
        df = df,
        years_in_era = NA_real_,
        from = label_a,
        to = label_b
      )
    }
  }
  list(
    effects = data.table::rbindlist(effects),
    comparisons = data.table::rbindlist(comparisons)
  )
}

policy_era_interaction <- function(counts = national_counts()) {
  groups <- list(
    list(label = "policy_2001_2006", selector = function(year) year <= 2006),
    list(label = "policy_2007_2019", selector = function(year) year >= 2007 && year <= 2019),
    list(label = "pandemic_recovery_2020_2024", selector = function(year) year >= 2020)
  )
  rows <- list()
  for (band_value in BANDS) {
    fitted <- fit_era_model(counts, band_value, groups)
    rows[[length(rows) + 1L]] <- fitted$effects
    rows[[length(rows) + 1L]] <- fitted$comparisons
  }
  data.table::rbindlist(rows, fill = TRUE)
}

period_interactions <- function(counts = national_counts()) {
  specifications <- list(
    list(
      model = "three_eras",
      groups = list(
        list(label = "pre_2001_2019", selector = function(year) year <= 2019),
        list(label = "covid_2020_2021", selector = function(year) year >= 2020 && year <= 2021),
        list(label = "recovery_2022_2024", selector = function(year) year >= 2022)
      )
    ),
    list(
      model = "policy_eras",
      groups = list(
        list(label = "policy_2001_2006", selector = function(year) year <= 2006),
        list(label = "policy_2007_2019", selector = function(year) year >= 2007 && year <= 2019),
        list(label = "pandemic_recovery_2020_2024", selector = function(year) year >= 2020)
      )
    )
  )
  rows <- list()
  for (specification in specifications) {
    for (band_value in BANDS) {
      fitted <- fit_era_model(counts, band_value, specification$groups)
      for (frame in list(fitted$effects, fitted$comparisons)) {
        exported <- data.table::copy(frame)
        exported[, model := specification$model]
        data.table::setcolorder(exported, c("model", setdiff(names(exported), "model")))
        rows[[length(rows) + 1L]] <- exported
      }
    }
  }
  data.table::rbindlist(rows, fill = TRUE)
}

arizona_negative_control_support <- function(counts = national_counts()) {
  rows <- list()
  for (endpoint_value in c(PEDESTRIAN, OTHER_FATAL, "all_fatal_crashes")) {
    for (band_value in BANDS) {
      subset <- counts[
        geo == "AZ_Phoenix" &
          width_minutes == 40 &
          endpoint == endpoint_value &
          band %chin% c(band_value, "night")
      ]
      table <- data.table::dcast(
        subset,
        year + season + sham ~ period + band,
        value.var = "events"
      )
      a <- table[[paste0("post_", band_value)]]
      b <- table$post_night
      c0 <- table[[paste0("pre_", band_value)]]
      d <- table$pre_night
      total <- a + b + c0 + d
      focal <- a + c0
      post <- a + b
      informative <- pmax(0, post - (total - focal)) < pmin(post, focal)
      row <- list(
        endpoint = endpoint_value,
        band = band_value,
        focal_events = sum(subset[band == band_value, events]),
        reference_events = sum(subset[band == "night", events]),
        candidate_tables = nrow(table),
        informative_tables = sum(informative),
        all_zero_tables = sum(total == 0)
      )
      for (specification in list(
        list(corrected = TRUE, label = "three_center_adjusted"),
        list(corrected = FALSE, label = "reduced_actual_date_only")
      )) {
        result <- try(
          estimate_allocation(
            counts,
            band = band_value,
            endpoint = endpoint_value,
            geography = "AZ_Phoenix",
            corrected = specification$corrected
          ),
          silent = TRUE
        )
        label <- specification$label
        if (inherits(result, "try-error")) {
          row[[paste0(label, "_status")]] <- paste0(
            "unsupported: ",
            sub("^Error[^:]*: ", "", as.character(result))
          )
          row[[paste0(label, "_ratio")]] <- NA_real_
          row[[paste0(label, "_low")]] <- NA_real_
          row[[paste0(label, "_high")]] <- NA_real_
        } else {
          row[[paste0(label, "_status")]] <- "fit"
          row[[paste0(label, "_ratio")]] <- result$ratio
          row[[paste0(label, "_low")]] <- result$low
          row[[paste0(label, "_high")]] <- result$high
        }
      }
      rows[[length(rows) + 1L]] <- data.table::as.data.table(row)
    }
  }
  data.table::rbindlist(rows, fill = TRUE)
}


alaska_inclusion_sensitivity <- function(counts = national_counts()) {
  alaska <- sensitivity_table("alaska_counts.csv.gz")
  keys <- c(
    "year", "season", "sham", "period", "width_minutes",
    "band", "endpoint", "window_days", "slot_minutes"
  )
  main <- counts[geo == "DST_main", c(keys, "events"), with = FALSE]
  combined <- data.table::rbindlist(list(main, alaska[, c(keys, "events"), with = FALSE]))[
    , .(events = sum(events)), by = keys
  ]
  combined[, geo := "DST_main_plus_Alaska"]
  data.table::setcolorder(combined, names(counts))
  expanded <- data.table::rbindlist(list(counts, combined), use.names = TRUE)

  rows <- list()
  for (geography_value in c("DST_main", "DST_main_plus_Alaska")) {
    for (band_value in BANDS) {
      fitted <- list()
      for (endpoint_value in c(PEDESTRIAN, OTHER_FATAL, "all_fatal_crashes")) {
        result <- estimate_allocation(
          expanded,
          band = band_value,
          endpoint = endpoint_value,
          geography = geography_value
        )
        fitted[[endpoint_value]] <- result
        rows[[length(rows) + 1L]] <- data.table::data.table(
          geography = geography_value,
          endpoint = endpoint_value,
          band = band_value,
          ratio = result$ratio,
          low = result$low,
          high = result$high,
          se = result$se,
          df = result$df,
          informative_strata = result$informative_strata,
          raw_strata = result$raw_strata
        )
      }
      direct <- contrast_allocation(
        fitted[[PEDESTRIAN]],
        fitted[[OTHER_FATAL]],
        paste0(geography_value, "_", band_value)
      )
      rows[[length(rows) + 1L]] <- data.table::data.table(
        geography = geography_value,
        endpoint = "direct",
        band = band_value,
        ratio = direct$ratio,
        low = direct$low,
        high = direct$high,
        se = direct$se,
        df = direct$df,
        informative_strata = NA_real_,
        raw_strata = NA_real_
      )
    }
  }
  data.table::rbindlist(rows)
}

review_diagnostics <- function(counts = national_counts()) {
  list(
    season_heterogeneity = seasonal_heterogeneity(counts),
    leave_one_year = leave_one_year_influence(counts),
    policy_era_interaction = policy_era_interaction(counts),
    period_interactions = period_interactions(counts),
    arizona_negative_control = arizona_negative_control_support(counts),
    alaska_inclusion = alaska_inclusion_sensitivity(counts)
  )
}
