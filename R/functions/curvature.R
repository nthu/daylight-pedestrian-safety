five_center_counts <- function() {
  base <- national_counts()
  early <- sensitivity_table("counts_placebo_-70.csv")
  late <- sensitivity_table("counts_placebo_+70.csv")
  keys <- c("year", "season", "period", "width_minutes", "band", "endpoint")

  main <- base[
    geo == "DST_main" & year %between% c(2002, 2023) & sham %in% c(-35, 0, 35),
    c("geo", keys, "sham", "events"),
    with = FALSE
  ]
  checks <- list(
    list(label = "early", frame = early, placebo = 35, main = -35),
    list(label = "late", frame = late, placebo = -35, main = 35)
  )
  for (check in checks) {
    x <- check$frame[sham == check$placebo, c(keys, "events"), with = FALSE]
    y <- main[sham == check$main, c(keys, "events"), with = FALSE]
    joined <- merge(x, y, by = keys, suffixes = c("_placebo", "_main"))
    if (nrow(joined) != nrow(x) || !identical(joined$events_placebo, joined$events_main)) {
      stop(check$label, " placebo calendar does not recover the shared main-calendar center.", call. = FALSE)
    }
  }

  outer_early <- early[sham == 0, c("geo", keys, "events"), with = FALSE]
  outer_early[, sham := -70L]
  outer_late <- late[sham == 0, c("geo", keys, "events"), with = FALSE]
  outer_late[, sham := 70L]
  columns <- c("geo", "year", "season", "sham", "period", "width_minutes", "band", "endpoint", "events")
  combined <- data.table::rbindlist(list(
    outer_early[, ..columns],
    main[, ..columns],
    outer_late[, ..columns]
  ))
  identity <- c("geo", "year", "season", "sham", "period", "width_minutes", "band", "endpoint")
  if (anyDuplicated(combined, by = identity)) stop("Duplicate five-center count identity.", call. = FALSE)
  if (!identical(sort(unique(combined$sham)), c(-70L, -35L, 0L, 35L, 70L))) {
    stop("Five-center support is incomplete.", call. = FALSE)
  }
  data.table::setorderv(combined, identity)
  combined
}

fit_five_center <- function(
    counts,
    band,
    endpoint,
    quadratic = TRUE,
    years = NULL) {
  band_value <- band
  endpoint_value <- endpoint
  selected <- counts[
    geo == "DST_main" &
      width_minutes == 40 &
      endpoint == endpoint_value &
      band %chin% c(band_value, "night")
  ]
  if (!is.null(years)) selected <- selected[year %in% years]
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
  zero_information <- table[!informative, .(year, season, sham)]
  table <- table[informative]
  strata <- nrow(table)
  if (strata < 10L) stop("Insufficient informative five-center strata.", call. = FALSE)

  nuisance <- list()
  extra_names <- character()
  for (year_value in sort(unique(table$year))) {
    for (season_value in c("autumn", "spring")) {
      nuisance[[length(nuisance) + 1L]] <- as.numeric(
        table$year == year_value & table$season == season_value
      )
      extra_names <- c(extra_names, paste0("post_", year_value, "_", season_value))
    }
  }
  center <- table$sham / 35
  for (season_value in c("autumn", "spring")) {
    mask <- as.numeric(table$season == season_value)
    nuisance[[length(nuisance) + 1L]] <- mask * center
    extra_names <- c(extra_names, paste0(season_value, "_linear"))
    if (quadratic) {
      nuisance[[length(nuisance) + 1L]] <- mask * center^2
      extra_names <- c(extra_names, paste0(season_value, "_quadratic"))
    }
  }
  direction <- ifelse(table$season == "spring", 1, -1) * ifelse(band == "evening", 1, -1)
  nuisance[[length(nuisance) + 1L]] <- as.numeric(table$sham == 0) * direction
  extra_names <- c(extra_names, "light_gain")
  extras <- do.call(cbind, nuisance)

  observations <- 2L * strata
  design <- matrix(0, nrow = observations, ncol = strata + ncol(extras))
  success <- failure <- numeric(observations)
  policy_year <- integer(observations)
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
    stop("Five-center grouped-binomial model failed to converge.", call. = FALSE)
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
    stop("Five-center information matrix is not full rank.", call. = FALSE)
  }
  bread <- solve(information)
  cluster_scores <- rowsum(score_rows, policy_year, reorder = TRUE)
  clusters <- as.integer(rownames(cluster_scores))
  G <- nrow(cluster_scores)
  N <- nrow(design)
  K <- ncol(design)
  correction <- G / (G - 1) * (N - 1) / (N - K)
  covariance <- bread %*% crossprod(cluster_scores) %*% bread * correction
  target <- match("light_gain", names)
  beta <- unname(fit$coefficients[target])
  standard_error <- sqrt(covariance[target, target])
  df <- G - 1L
  interval <- t_interval(beta, standard_error, df)
  influence <- as.vector(cluster_scores %*% bread[, target]) * sqrt(correction)
  names(influence) <- clusters
  if (!isTRUE(all.equal(sum(influence^2), standard_error^2, tolerance = 1e-8))) {
    stop("Five-center policy-year influence decomposition mismatch.", call. = FALSE)
  }

  quadratic_indices <- if (quadratic) {
    match(paste0(c("autumn", "spring"), "_quadratic"), names)
  } else integer()
  if (length(quadratic_indices)) {
    coefficients <- fit$coefficients[quadratic_indices]
    variance <- covariance[quadratic_indices, quadratic_indices, drop = FALSE]
    wald <- as.numeric(crossprod(coefficients, solve(variance, coefficients)))
    quadratic_p <- stats::pf(wald / length(quadratic_indices), length(quadratic_indices), df, lower.tail = FALSE)
  } else {
    wald <- quadratic_p <- NA_real_
  }

  c(
    list(
      band = band,
      endpoint = endpoint,
      model = if (quadratic) "five_center_quadratic" else "five_center_linear"
    ),
    interval,
    list(
      clusters = G,
      allocation_observations = N,
      parameters = K,
      informative_strata = strata,
      zero_information_strata = nrow(zero_information),
      cr1_factor = correction,
      information_condition = kappa(information),
      quadratic_wald = wald,
      quadratic_df_num = if (quadratic) 2L else NA_integer_,
      quadratic_df_denom = if (quadratic) df else NA_integer_,
      quadratic_p_value = quadratic_p,
      autumn_quadratic = if (quadratic) fit$coefficients[match("autumn_quadratic", names)] else NA_real_,
      spring_quadratic = if (quadratic) fit$coefficients[match("spring_quadratic", names)] else NA_real_,
      year_influence = influence
    )
  )
}

five_center_analysis <- function(full = TRUE) {
  counts <- five_center_counts()
  estimates <- list()
  deletions <- list()
  fitted <- list()

  for (quadratic in c(FALSE, TRUE)) {
    model <- if (quadratic) "five_center_quadratic" else "five_center_linear"
    for (endpoint in ENDPOINTS) {
      for (band in BANDS) {
        result <- fit_five_center(counts, band, endpoint, quadratic = quadratic)
        fitted[[paste(model, endpoint, band, sep = "|")]] <- result
        row <- result[setdiff(names(result), "year_influence")]
        estimates[[length(estimates) + 1L]] <- data.table::as.data.table(row)
      }
    }
    for (band in BANDS) {
      pedestrian <- fitted[[paste(model, PEDESTRIAN, band, sep = "|")]]
      other <- fitted[[paste(model, OTHER_FATAL, band, sep = "|")]]
      years <- intersect(names(pedestrian$year_influence), names(other$year_influence))
      influence <- pedestrian$year_influence[years] - other$year_influence[years]
      beta <- pedestrian$log_ratio - other$log_ratio
      se <- sqrt(sum(influence^2))
      interval <- t_interval(beta, se, length(years) - 1L)
      estimates[[length(estimates) + 1L]] <- data.table::data.table(
        band = band,
        endpoint = "direct",
        model = model,
        log_ratio = interval$log_ratio,
        ratio = interval$ratio,
        se = interval$se,
        low = interval$low,
        high = interval$high,
        p_value = interval$p_value,
        df = interval$df,
        clusters = length(years),
        allocation_observations = NA_real_,
        parameters = NA_real_,
        informative_strata = NA_real_,
        zero_information_strata = NA_real_,
        cr1_factor = NA_real_,
        information_condition = NA_real_,
        quadratic_wald = NA_real_,
        quadratic_df_num = NA_real_,
        quadratic_df_denom = NA_real_,
        quadratic_p_value = NA_real_,
        autumn_quadratic = NA_real_,
        spring_quadratic = NA_real_
      )
    }
  }

  if (full) {
    for (endpoint in ENDPOINTS) {
      for (band in BANDS) {
        full_result <- fitted[[paste("five_center_quadratic", endpoint, band, sep = "|")]]
        deleted <- numeric()
        for (omitted in 2002:2023) {
          result <- fit_five_center(
            counts,
            band,
            endpoint,
            quadratic = TRUE,
            years = setdiff(2002:2023, omitted)
          )
          deleted <- c(deleted, result$log_ratio)
          deletions[[length(deletions) + 1L]] <- data.table::data.table(
            band = band,
            endpoint = endpoint,
            omitted_policy_year = omitted,
            log_ratio = result$log_ratio,
            ratio = result$ratio
          )
        }
        interval <- jackknife_interval(full_result$log_ratio, deleted)
        estimates[[length(estimates) + 1L]] <- data.table::data.table(
          band = band,
          endpoint = endpoint,
          model = "five_center_quadratic_year_jackknife",
          log_ratio = interval$log_ratio,
          ratio = interval$ratio,
          se = interval$se,
          low = interval$low,
          high = interval$high,
          p_value = interval$p_value,
          df = interval$df,
          clusters = interval$blocks
        )
      }
    }

    deleted_table <- data.table::rbindlist(deletions)
    for (band in BANDS) {
      band_value <- band
      pedestrian <- fitted[[paste("five_center_quadratic", PEDESTRIAN, band, sep = "|")]]
      other <- fitted[[paste("five_center_quadratic", OTHER_FATAL, band, sep = "|")]]
      ped <- deleted_table[band == band_value & endpoint == PEDESTRIAN][order(omitted_policy_year)]
      oth <- deleted_table[band == band_value & endpoint == OTHER_FATAL][order(omitted_policy_year)]
      values <- ped$log_ratio - oth$log_ratio
      interval <- jackknife_interval(pedestrian$log_ratio - other$log_ratio, values)
      estimates[[length(estimates) + 1L]] <- data.table::data.table(
        band = band,
        endpoint = "direct",
        model = "five_center_quadratic_year_jackknife",
        log_ratio = interval$log_ratio,
        ratio = interval$ratio,
        se = interval$se,
        low = interval$low,
        high = interval$high,
        p_value = interval$p_value,
        df = interval$df,
        clusters = interval$blocks
      )
    }
  }

  list(
    five_center_curvature = data.table::rbindlist(estimates, fill = TRUE),
    five_center_deletions = if (length(deletions)) data.table::rbindlist(deletions) else data.table::data.table()
  )
}

exact_rademacher_multiplier <- function(log_ratio, influences, chunk_size = 2^17) {
  values <- as.numeric(influences)
  clusters <- length(values)
  if (clusters < 3L || any(!is.finite(values))) {
    stop("At least three finite policy-year influences are required.", call. = FALSE)
  }
  if (abs(sum(values)) > 1e-8) {
    stop("Influence contributions must sum to zero at the fitted model.", call. = FALSE)
  }
  variance <- sum(values^2)
  standard_error <- sqrt(variance)
  observed_t <- abs(log_ratio / standard_error)
  enumerated <- 2^(clusters - 1L)
  statistics <- numeric(enumerated)
  position <- 1L
  for (start in seq.int(0L, enumerated - 1L, by = chunk_size)) {
    finish <- min(start + chunk_size - 1L, enumerated - 1L)
    identifiers <- as.integer(start:finish)
    signs <- matrix(1, nrow = length(identifiers), ncol = clusters)
    for (column in 2:clusters) {
      signs[, column] <- 1 - 2 * bitwAnd(bitwShiftR(identifiers, column - 2L), 1L)
    }
    sums <- as.vector(signs %*% values)
    denominator <- sqrt(pmax(variance - sums^2 / clusters, 1e-14))
    current <- abs(sums / denominator)
    statistics[position:(position + length(current) - 1L)] <- current
    position <- position + length(current)
  }
  critical <- unname(stats::quantile(statistics, 0.95, type = 1))
  p_value <- mean(statistics >= observed_t - 1e-12)
  list(
    wild_p_value = p_value,
    wild_critical_95 = critical,
    wild_low = exp(log_ratio - critical * standard_error),
    wild_high = exp(log_ratio + critical * standard_error),
    observed_t = observed_t,
    clusters = clusters,
    unique_sign_patterns = enumerated,
    equivalent_full_sign_patterns = 2^clusters
  )
}

wild_cluster_calibration <- function(counts = national_counts()) {
  fitted <- list()
  rows <- list()
  for (endpoint in ENDPOINTS) {
    for (band in BANDS) {
      fitted[[paste(endpoint, band, sep = "|")]] <- estimate_allocation(
        counts,
        band = band,
        endpoint = endpoint
      )
    }
  }
  for (band in BANDS) {
    pedestrian <- fitted[[paste(PEDESTRIAN, band, sep = "|")]]
    calibrated <- exact_rademacher_multiplier(
      pedestrian$log_ratio,
      pedestrian$year_influence
    )
    rows[[length(rows) + 1L]] <- data.table::as.data.table(c(
      list(
        comparison = paste0("pedestrian_", band),
        band = band,
        endpoint = PEDESTRIAN,
        ratio = pedestrian$ratio,
        cr1_low = pedestrian$low,
        cr1_high = pedestrian$high,
        cr1_p_value = pedestrian$p_value
      ),
      calibrated
    ))

    other <- fitted[[paste(OTHER_FATAL, band, sep = "|")]]
    years <- intersect(names(pedestrian$year_influence), names(other$year_influence))
    influence <- pedestrian$year_influence[years] - other$year_influence[years]
    beta <- pedestrian$log_ratio - other$log_ratio
    se <- sqrt(sum(influence^2))
    df <- length(years) - 1L
    critical <- stats::qt(0.975, df)
    calibrated <- exact_rademacher_multiplier(beta, influence)
    rows[[length(rows) + 1L]] <- data.table::as.data.table(c(
      list(
        comparison = paste0("pedestrian_vs_other_", band),
        band = band,
        endpoint = "direct",
        ratio = exp(beta),
        cr1_low = exp(beta - critical * se),
        cr1_high = exp(beta + critical * se),
        cr1_p_value = 2 * stats::pt(-abs(beta / se), df)
      ),
      calibrated
    ))
  }
  data.table::rbindlist(rows, fill = TRUE)
}
