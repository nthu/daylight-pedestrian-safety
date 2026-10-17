build_allocation_frame <- function(
    counts,
    band = "evening",
    reference = "night",
    endpoint = "fatal_pedestrian_crashes",
    geography = "DST_main",
    width_minutes = 40,
    years = NULL,
    seasons = NULL,
    corrected = TRUE) {
  band_value <- band
  reference_value <- reference
  endpoint_value <- endpoint
  geography_value <- geography
  width_value <- width_minutes
  selected <- counts[
    geo == geography_value &
      width_minutes == width_value &
      endpoint == endpoint_value &
      band %chin% c(band_value, reference_value)
  ]
  if (!is.null(years)) selected <- selected[year %in% years]
  if (!is.null(seasons)) selected <- selected[season %chin% seasons]
  if (!corrected) selected <- selected[sham == 0]

  table <- data.table::dcast(
    selected,
    year + season + sham ~ period + band,
    value.var = "events"
  )
  data.table::setorderv(table, c("year", "season", "sham"))

  required <- c(
    paste0("post_", band),
    paste0("post_", reference),
    paste0("pre_", band),
    paste0("pre_", reference)
  )
  if (!all(required %in% names(table))) {
    stop("The focal/reference table is incomplete.", call. = FALSE)
  }

  a <- table[[paste0("post_", band)]]
  b <- table[[paste0("post_", reference)]]
  c0 <- table[[paste0("pre_", band)]]
  d <- table[[paste0("pre_", reference)]]
  total <- a + b + c0 + d
  focal <- a + c0
  post <- a + b
  informative <- pmax(0, post - (total - focal)) < pmin(post, focal)
  zero_information <- table[!informative, .(year, season, sham)]
  table <- table[informative]

  strata <- nrow(table)
  if (strata < 6L) stop("Insufficient informative strata.", call. = FALSE)

  nuisance <- list()
  extra_names <- character()
  available_seasons <- sort(unique(table$season))
  if (corrected) {
    for (year_value in sort(unique(table$year))) {
      for (season_value in available_seasons) {
        nuisance[[length(nuisance) + 1L]] <- as.numeric(
          table$year == year_value & table$season == season_value
        )
        extra_names <- c(extra_names, paste0("post_", year_value, "_", season_value))
      }
    }
    for (season_value in available_seasons) {
      nuisance[[length(nuisance) + 1L]] <-
        as.numeric(table$season == season_value) * table$sham / 35
      extra_names <- c(extra_names, paste0("post_", season_value, "_centre_slope"))
    }
  }

  direction <- ifelse(table$season == "spring", 1, -1) *
    ifelse(startsWith(band, "evening"), 1, -1)
  nuisance[[length(nuisance) + 1L]] <- as.numeric(table$sham == 0) * direction
  extra_names <- c(extra_names, "light_gain")
  extras <- do.call(cbind, nuisance)

  rows <- 2L * strata
  design <- matrix(0, nrow = rows, ncol = strata + ncol(extras))
  success <- failure <- numeric(rows)
  policy_year <- integer(rows)
  row <- 1L
  for (index in seq_len(strata)) {
    for (post_period in 0:1) {
      design[row, index] <- 1
      if (post_period == 1L) {
        design[row, (strata + 1L):ncol(design)] <- extras[index, ]
      }
      period <- if (post_period == 0L) "pre" else "post"
      success[row] <- table[[paste0(period, "_", band)]][index]
      failure[row] <- table[[paste0(period, "_", reference)]][index]
      policy_year[row] <- table$year[index]
      row <- row + 1L
    }
  }

  column_names <- c(paste0("stratum_", seq_len(strata) - 1L), extra_names)
  colnames(design) <- column_names
  list(
    table = table,
    design = design,
    success = success,
    failure = failure,
    policy_year = policy_year,
    names = column_names,
    zero_information = zero_information
  )
}

estimate_allocation <- function(
    counts,
    band = "evening",
    reference = "night",
    endpoint = "fatal_pedestrian_crashes",
    geography = "DST_main",
    width_minutes = 40,
    years = NULL,
    seasons = NULL,
    corrected = TRUE) {
  built <- build_allocation_frame(
    counts,
    band = band,
    reference = reference,
    endpoint = endpoint,
    geography = geography,
    width_minutes = width_minutes,
    years = years,
    seasons = seasons,
    corrected = corrected
  )

  fit <- stats::glm.fit(
    x = built$design,
    y = cbind(built$success, built$failure),
    family = stats::binomial(),
    control = stats::glm.control(epsilon = 1e-10, maxit = 200)
  )
  if (!isTRUE(fit$converged) || any(!is.finite(fit$coefficients))) {
    stop("Allocation model did not converge to finite coefficients.", call. = FALSE)
  }

  total <- built$success + built$failure
  probability <- fit$fitted.values
  residual <- built$success - total * probability
  score_rows <- built$design * residual
  information <- crossprod(
    built$design,
    built$design * (total * probability * (1 - probability))
  )
  eigenvalues <- eigen((information + t(information)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) <= max(eigenvalues) * 1e-12) {
    stop("Allocation information matrix is not full rank under the fixed tolerance.", call. = FALSE)
  }
  bread <- solve(information)
  cluster_scores <- rowsum(score_rows, built$policy_year, reorder = TRUE)

  clusters <- nrow(cluster_scores)
  observations <- nrow(built$design)
  parameters <- ncol(built$design)
  if (clusters < 2L || observations <= parameters) {
    stop("Insufficient clusters or residual allocation observations.", call. = FALSE)
  }
  correction <- clusters / (clusters - 1) *
    (observations - 1) / (observations - parameters)
  covariance <- bread %*% crossprod(cluster_scores) %*% bread * correction

  target <- match("light_gain", built$names)
  log_ratio <- unname(fit$coefficients[target])
  standard_error <- sqrt(covariance[target, target])
  interval <- t_interval(log_ratio, standard_error, clusters - 1L)
  influence <- as.vector(cluster_scores %*% bread[, target]) * sqrt(correction)
  names(influence) <- rownames(cluster_scores)
  if (!isTRUE(all.equal(sum(influence^2), standard_error^2, tolerance = 1e-8))) {
    stop("Policy-year influence decomposition does not match the variance.", call. = FALSE)
  }

  c(
    interval,
    list(
      band = band,
      reference = reference,
      endpoint = endpoint,
      geography = geography,
      width_minutes = width_minutes,
      corrected = corrected,
      raw_strata = nrow(built$table) + nrow(built$zero_information),
      informative_strata = nrow(built$table),
      zero_information_strata = built$zero_information,
      years = sort(unique(built$policy_year)),
      year_influence = influence,
      year_CR1_factor = correction,
      coefficients = stats::setNames(unname(fit$coefficients), built$names),
      covariance = covariance,
      design_columns = built$names,
      input_events_across_reused_tables = sum(built$success + built$failure),
      converged = fit$converged
    )
  )
}

contrast_allocation <- function(left, right, name) {
  shared <- intersect(names(left$year_influence), names(right$year_influence))
  cross_covariance <- sum(
    left$year_influence[shared] * right$year_influence[shared]
  )
  variance <- left$se^2 + right$se^2 - 2 * cross_covariance
  if (variance < -1e-10) stop("Negative direct-comparison variance.", call. = FALSE)
  standard_error <- sqrt(max(variance, 0))
  log_ratio <- left$log_ratio - right$log_ratio
  degrees_freedom <- min(left$df, right$df)
  c(
    list(contrast = name),
    t_interval(log_ratio, standard_error, degrees_freedom),
    list(
      shared_years = length(shared),
      cross_covariance = unname(cross_covariance)
    )
  )
}
