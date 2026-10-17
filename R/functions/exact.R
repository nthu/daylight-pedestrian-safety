log_sum_exp <- function(x) {
  maximum <- max(x)
  maximum + log(sum(exp(x - maximum)))
}

prepare_exact_tables <- function(cells, design, years, names) {
  cells <- as.matrix(cells)
  storage.mode(cells) <- "integer"
  design <- as.matrix(design)
  if (ncol(cells) != 4L || nrow(cells) != nrow(design)) {
    stop("Expected n-by-4 tables and a matching design matrix.", call. = FALSE)
  }
  if (any(cells < 0) || any(!is.finite(cells)) || any(!is.finite(design))) {
    stop("Counts and design values must be finite and nonnegative.", call. = FALSE)
  }
  a <- cells[, 1]
  post_total <- cells[, 1] + cells[, 2]
  pre_total <- cells[, 3] + cells[, 4]
  focal_total <- cells[, 1] + cells[, 3]
  lower <- pmax(0L, focal_total - pre_total)
  upper <- pmin(post_total, focal_total)
  support <- Map(seq.int, lower, upper)
  log_weight <- Map(
    function(x, p, q, f) lchoose(p, x) + lchoose(q, f - x),
    support,
    post_total,
    pre_total,
    focal_total
  )
  list(
    cells = cells,
    design = design,
    years = as.integer(years),
    names = names,
    a = a,
    post_total = post_total,
    pre_total = pre_total,
    focal_total = focal_total,
    lower = lower,
    upper = upper,
    informative = upper > lower,
    support = support,
    log_weight = log_weight,
    observed_log_weight = lchoose(post_total, a) + lchoose(pre_total, focal_total - a)
  )
}

exact_evaluate <- function(model, beta) {
  eta <- as.vector(model$design %*% beta)
  n <- length(eta)
  means <- variances <- residuals <- numeric(n)
  log_likelihood <- 0
  for (index in seq_len(n)) {
    values <- model$support[[index]]
    log_terms <- model$log_weight[[index]] + values * eta[index]
    normalizer <- log_sum_exp(log_terms)
    probabilities <- exp(log_terms - normalizer)
    means[index] <- sum(probabilities * values)
    variances[index] <- sum(probabilities * (values - means[index])^2)
    residuals[index] <- model$a[index] - means[index]
    log_likelihood <- log_likelihood +
      model$observed_log_weight[index] + model$a[index] * eta[index] - normalizer
  }
  scores <- model$design * residuals
  list(
    log_likelihood = log_likelihood,
    score = colSums(scores),
    information = crossprod(model$design, model$design * variances),
    scores = scores,
    means = means,
    variances = variances
  )
}

fit_exact_tables <- function(model, start = NULL, tolerance = 2e-8, max_iterations = 100L) {
  beta <- if (is.null(start)) numeric(ncol(model$design)) else as.numeric(start)
  converged <- FALSE
  for (iteration in seq_len(max_iterations)) {
    evaluated <- exact_evaluate(model, beta)
    error <- max(abs(evaluated$score))
    if (error < tolerance) {
      converged <- TRUE
      break
    }
    eigenvalues <- eigen(evaluated$information, symmetric = TRUE, only.values = TRUE)$values
    if (min(eigenvalues) <= max(eigenvalues) * 1e-12) {
      stop("Exact-likelihood information matrix is not full rank.", call. = FALSE)
    }
    step <- solve(evaluated$information, evaluated$score)
    scale <- 1
    accepted <- FALSE
    for (attempt in seq_len(35L)) {
      candidate <- beta + scale * step
      candidate_value <- exact_evaluate(model, candidate)
      if (candidate_value$log_likelihood >= evaluated$log_likelihood - 1e-10 &&
          max(abs(candidate_value$score)) < error) {
        beta <- candidate
        accepted <- TRUE
        break
      }
      scale <- scale / 2
    }
    if (!accepted) stop("Exact-likelihood line search failed.", call. = FALSE)
  }
  if (!converged) stop("Exact-likelihood iteration limit reached.", call. = FALSE)

  evaluated <- exact_evaluate(model, beta)
  information <- evaluated$information
  bread <- solve(information)
  clusters <- sort(unique(model$years))
  cluster_scores <- do.call(
    rbind,
    lapply(clusters, function(year) colSums(evaluated$scores[model$years == year, , drop = FALSE]))
  )
  covariance <- bread %*% crossprod(cluster_scores) %*% bread *
    length(clusters) / (length(clusters) - 1)
  target <- match("light_gain", model$names)
  standard_error <- sqrt(covariance[target, target])
  interval <- t_interval(beta[target], standard_error, length(clusters) - 1L)

  c(
    interval,
    list(
      beta = stats::setNames(beta, model$names),
      covariance = covariance,
      working_information_se = sqrt(bread[target, target]),
      years = clusters,
      tables_total = nrow(model$cells),
      tables_informative = sum(model$informative),
      tables_degenerate = sum(!model$informative),
      score_max = max(abs(evaluated$score)),
      information_condition = kappa(information),
      log_likelihood = evaluated$log_likelihood,
      iterations = iteration
    )
  )
}

make_exact_model <- function(frame, band, years = NULL, seasons = NULL) {
  band_value <- band
  selected <- data.table::copy(frame[band == band_value])
  if (!is.null(years)) selected <- selected[year %in% years]
  if (!is.null(seasons)) selected <- selected[season %chin% seasons]
  if (anyDuplicated(selected[, .(year, season, sham)])) {
    stop("Duplicated exact-table identity.", call. = FALSE)
  }
  data.table::setorderv(selected, c("year", "season", "sham"))
  triplets <- selected[, .(centers = paste(sort(sham), collapse = ",")), by = .(year, season)]
  if (any(triplets$centers != "-35,0,35")) {
    stop("Every year-season requires real and both sham centres.", call. = FALSE)
  }

  year_seasons <- unique(selected[, .(year, season)])
  names <- paste0("yearseason_", year_seasons$year, "_", year_seasons$season)
  columns <- lapply(
    seq_len(nrow(year_seasons)),
    function(index) as.numeric(
      selected$year == year_seasons$year[index] &
        selected$season == year_seasons$season[index]
    )
  )
  for (season_value in sort(unique(selected$season))) {
    names <- c(names, paste0("slope_", season_value))
    columns[[length(columns) + 1L]] <-
      as.numeric(selected$season == season_value) * selected$sham / 35
  }
  direction <- ifelse(selected$season == "spring", 1, -1) *
    ifelse(startsWith(band, "evening"), 1, -1)
  names <- c(names, "light_gain")
  columns[[length(columns) + 1L]] <- as.numeric(selected$sham == 0) * direction

  cells <- as.matrix(
    selected[, .(focal_post, dark_post, focal_pre, dark_pre)]
  )
  model <- prepare_exact_tables(
    cells,
    do.call(cbind, columns),
    selected$year,
    names
  )
  list(model = model, frame = selected)
}

exact_block_intervals <- function(frame, band, block_sizes = 1:4) {
  built <- make_exact_model(frame, band)
  full <- fit_exact_tables(built$model)
  years <- 2001:2024
  rows <- list()
  for (block_size in block_sizes) {
    deleted <- numeric()
    for (start in seq(2001, 2024, by = block_size)) {
      omitted <- start:(start + block_size - 1L)
      retained <- setdiff(years, omitted)
      reduced <- make_exact_model(frame, band, years = retained)
      deleted <- c(deleted, fit_exact_tables(reduced$model)$log_ratio)
    }
    interval <- jackknife_interval(full$log_ratio, deleted)
    rows[[length(rows) + 1L]] <- data.table::data.table(
      band = band,
      block_size_years = block_size,
      blocks = interval$blocks,
      df = interval$df,
      log_ratio = interval$log_ratio,
      ratio = interval$ratio,
      se = interval$se,
      low = interval$low,
      high = interval$high,
      p_value = interval$p_value
    )
  }
  data.table::rbindlist(rows)
}
