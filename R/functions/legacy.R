estimate_legacy_weeks <- function() {
  counts <- sensitivity_table("legacy_weeks.csv")
  if (nrow(counts) != 1152L || any(counts$events < 0)) {
    stop("Unexpected legacy-week input.", call. = FALSE)
  }
  estimates <- list()
  pretests <- list()

  for (band in BANDS) {
    band_value <- band
    selected <- counts[band %chin% c(band_value, "night")]
    table <- data.table::dcast(
      selected,
      year + season + week ~ band,
      value.var = "events"
    )
    data.table::setorderv(table, c("year", "season", "week"))
    unit <- interaction(table$year, table$season, drop = TRUE)
    unit_levels <- levels(unit)
    design <- stats::model.matrix(~ 0 + unit)
    names <- paste0("baseline_", unit_levels)

    extra <- list()
    for (season_value in SEASONS) {
      for (week_value in c(-4, -3, -2, 0, 1, 2, 3)) {
        extra[[length(extra) + 1L]] <- as.numeric(
          table$season == season_value & table$week == week_value
        )
        names <- c(names, paste0(season_value, "_week_", week_value))
      }
    }
    design <- cbind(design, do.call(cbind, extra))
    colnames(design) <- names
    success <- table[[band]]
    failure <- table$night

    fit <- stats::glm.fit(
      x = design,
      y = cbind(success, failure),
      family = stats::binomial(),
      control = stats::glm.control(epsilon = 1e-10, maxit = 200)
    )
    if (!fit$converged || any(!is.finite(fit$coefficients))) {
      stop("Legacy event-week model did not converge.", call. = FALSE)
    }
    total <- success + failure
    residual <- success - total * fit$fitted.values
    score_rows <- design * residual
    information <- crossprod(
      design,
      design * (total * fit$fitted.values * (1 - fit$fitted.values))
    )
    bread <- solve(information)
    cluster_scores <- rowsum(score_rows, table$year, reorder = TRUE)
    clusters <- nrow(cluster_scores)
    used_groups <- sum(total > 0)
    correction <- clusters / (clusters - 1) *
      (used_groups - 1) / (used_groups - ncol(design))
    covariance <- bread %*% crossprod(cluster_scores) %*% bread * correction
    critical <- stats::qt(0.975, clusters - 1)

    for (season_value in SEASONS) {
      for (week_value in c(-4, -3, -2, 0, 1, 2, 3)) {
        term <- paste0(season_value, "_week_", week_value)
        index <- match(term, names)
        beta <- fit$coefficients[index]
        standard_error <- sqrt(covariance[index, index])
        estimates[[length(estimates) + 1L]] <- data.table::data.table(
          band = band,
          term = term,
          season = season_value,
          week = week_value,
          log_ratio = beta,
          ratio = exp(beta),
          se = standard_error,
          low = exp(beta - critical * standard_error),
          high = exp(beta + critical * standard_error),
          p_value = 2 * stats::pt(-abs(beta / standard_error), clusters - 1),
          df = clusters - 1
        )
      }
      indices <- match(
        paste0(season_value, "_week_", c(-4, -3, -2)),
        names
      )
      beta <- fit$coefficients[indices]
      covariance_pre <- covariance[indices, indices, drop = FALSE]
      if (qr(covariance_pre)$rank != 3L) {
        stop("Pre-period covariance is not full rank.", call. = FALSE)
      }
      statistic <- as.numeric(crossprod(beta, solve(covariance_pre, beta)) / 3)
      pretests[[length(pretests) + 1L]] <- data.table::data.table(
        band = band,
        season = season_value,
        F = statistic,
        df_num = 3L,
        df_denom = clusters - 1L,
        p_value = stats::pf(statistic, 3, clusters - 1, lower.tail = FALSE)
      )
    }
  }
  list(
    estimates = data.table::rbindlist(estimates),
    pretests = data.table::rbindlist(pretests)
  )
}
