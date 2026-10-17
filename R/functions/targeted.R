make_triplet <- function(cells, year = 2001L, season = "spring", direction = 1) {
  cells <- as.matrix(cells)
  storage.mode(cells) <- "integer"
  if (!identical(dim(cells), c(3L, 4L)) || any(cells < 0)) {
    stop("Three nonnegative integer 2-by-2 tables are required.", call. = FALSE)
  }
  post_total <- cells[, 1] + cells[, 2]
  pre_total <- cells[, 3] + cells[, 4]
  focal_total <- cells[, 1] + cells[, 3]
  lower <- pmax(0L, focal_total - pre_total)
  upper <- pmin(post_total, focal_total)
  total_focal_post <- sum(cells[, 1])
  end_difference <- cells[3, 1] - cells[1, 1]

  support <- numeric()
  log_weight <- numeric()
  for (middle in seq.int(lower[2], upper[2])) {
    if ((total_focal_post - middle - end_difference) %% 2 != 0) next
    early <- (total_focal_post - middle - end_difference) / 2
    late <- (total_focal_post - middle + end_difference) / 2
    values <- c(early, middle, late)
    if (any(values < lower) || any(values > upper)) next
    support <- c(support, middle)
    log_weight <- c(
      log_weight,
      sum(lchoose(post_total, values) + lchoose(pre_total, focal_total - values))
    )
  }
  if (!length(support)) stop("Conditioned triplet has empty support.", call. = FALSE)
  list(
    year = as.integer(year),
    season = season,
    direction = as.numeric(direction),
    observed = unname(cells[2, 1]),
    support = support,
    log_weight = log_weight,
    T = total_focal_post,
    U = end_difference
  )
}

triplet_moments <- function(triplet, beta) {
  log_terms <- triplet$log_weight + triplet$direction * beta * triplet$support
  normalizer <- log_sum_exp(log_terms)
  probability <- exp(log_terms - normalizer)
  mean <- sum(probability * triplet$support)
  variance <- sum(probability * (triplet$support - mean)^2)
  observed_index <- match(triplet$observed, triplet$support)
  if (is.na(observed_index)) stop("Observed triplet count is absent from support.", call. = FALSE)
  c(
    log_likelihood = log_terms[observed_index] - normalizer,
    score = unname(triplet$direction * (triplet$observed - mean)),
    information = unname(triplet$direction^2 * variance)
  )
}

build_triplets <- function(counts, band, endpoint, years = NULL) {
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
  groups <- split(table, interaction(table$year, table$season, drop = TRUE))
  lapply(groups, function(group) {
    group <- group[match(c(-35, 0, 35), group$sham)]
    if (anyNA(group$sham)) stop("Incomplete year-season triplet.", call. = FALSE)
    cells <- as.matrix(group[, .SD, .SDcols = c(
      paste0("post_", band),
      "post_night",
      paste0("pre_", band),
      "pre_night"
    )])
    direction <- ifelse(group$season[1] == "spring", 1, -1) *
      ifelse(band == "evening", 1, -1)
    make_triplet(cells, group$year[1], group$season[1], direction)
  })
}

fit_triplets <- function(triplets) {
  informative <- triplets[vapply(triplets, function(x) length(x$support) > 1L, logical(1))]
  if (length(informative) < 3L) stop("Too few informative triplets.", call. = FALSE)
  sum_moments <- function(beta) {
    Reduce(`+`, lapply(informative, triplet_moments, beta = beta))
  }
  left <- sum_moments(-40)["score"]
  right <- sum_moments(40)["score"]
  if (!(left > 0 && right < 0)) {
    stop("Finite triplet effect is not supported.", call. = FALSE)
  }
  beta <- stats::uniroot(
    function(value) sum_moments(value)["score"],
    lower = -40,
    upper = 40,
    tol = 1e-12
  )$root
  moments <- sum_moments(beta)
  if (abs(moments["score"]) > 1e-7 || moments["information"] <= 0) {
    stop("Triplet score or information check failed.", call. = FALSE)
  }
  list(
    log_ratio = unname(beta),
    ratio = unname(exp(beta)),
    score = unname(moments["score"]),
    information = unname(moments["information"]),
    informative_triplets = length(informative),
    total_triplets = length(triplets),
    degenerate_triplets = length(triplets) - length(informative),
    log_likelihood = unname(moments["log_likelihood"])
  )
}

standardized_inputs <- function(counts, band) {
  band_value <- band
  selected <- counts[
    geo == "DST_main" &
      width_minutes == 40 &
      endpoint %chin% ENDPOINTS &
      band %chin% c(band_value, "night")
  ]
  table <- data.table::dcast(
    selected,
    year + season + endpoint + sham + period ~ band,
    value.var = "events"
  )
  table[, total := get(band_value) + night]
  table[, fraction := get(band_value) / total]
  support <- table[, .(complete = .N == 24L && all(total > 0)), by = year]
  years <- support[complete == TRUE, year]
  excluded <- support[complete == FALSE, year]
  table <- table[year %in% years]
  table[, direction := ifelse(season == "spring", 1, -1) *
    ifelse(band_value == "evening", 1, -1)]
  table[, gain_side := (period == "post" & direction == 1) |
    (period == "pre" & direction == -1)]
  list(frame = table, years = years, excluded = excluded)
}

standardized_estimate <- function(frame, years) {
  selected <- frame[year %in% years]
  expected_size <- 2L * length(years)
  average <- selected[, .(fraction = mean(fraction), size = .N),
    by = .(endpoint, sham, gain_side)]
  if (any(average$size != expected_size)) {
    stop("Common temporal support is not equal.", call. = FALSE)
  }
  contrast <- list()
  parts <- list()
  for (endpoint_value in ENDPOINTS) {
    center_contrast <- numeric()
    for (center in c(-35, 0, 35)) {
      gain <- average[endpoint == endpoint_value & sham == center & gain_side == TRUE, fraction]
      loss <- average[endpoint == endpoint_value & sham == center & gain_side == FALSE, fraction]
      if (length(gain) != 1L || length(loss) != 1L ||
          gain <= 0 || gain >= 1 || loss <= 0 || loss >= 1) {
        stop("Boundary or incomplete marginal probability.", call. = FALSE)
      }
      value <- stats::qlogis(gain) - stats::qlogis(loss)
      center_contrast[as.character(center)] <- value
      coefficient <- if (center == 0) 1 else -0.5
      parts[[length(parts) + 1L]] <- data.table::data.table(
        endpoint = endpoint_value,
        sham = center,
        gain_fraction = gain,
        loss_fraction = loss,
        log_contrast = value,
        contrast = exp(value),
        coefficient_in_symmetric_sum = coefficient,
        contribution = coefficient * value,
        weight_per_year_season = 1 / expected_size
      )
    }
    contrast[[endpoint_value]] <- center_contrast["0"] -
      0.5 * (center_contrast["-35"] + center_contrast["35"])
  }
  contrast$direct <- contrast[[PEDESTRIAN]] - contrast[[OTHER_FATAL]]
  list(estimate = contrast, parts = data.table::rbindlist(parts))
}

triplet_analysis <- function(counts) {
  years <- 2001:2024
  estimates <- list()
  support <- list()
  stored <- list()
  for (band in BANDS) {
    for (endpoint in ENDPOINTS) {
      triplets <- build_triplets(counts, band, endpoint)
      full <- fit_triplets(triplets)
      deleted <- vapply(
        years,
        function(year) fit_triplets(
          triplets[vapply(triplets, function(x) x$year != year, logical(1))]
        )$log_ratio,
        numeric(1)
      )
      interval <- jackknife_interval(full$log_ratio, deleted)
      estimates[[length(estimates) + 1L]] <- data.table::data.table(
        band = band,
        endpoint = endpoint,
        ratio = interval$ratio,
        log_ratio = interval$log_ratio,
        se = interval$se,
        low = interval$low,
        high = interval$high,
        p_value = interval$p_value,
        df = interval$df,
        informative_triplets = full$informative_triplets,
        total_triplets = full$total_triplets,
        degenerate_triplets = full$degenerate_triplets
      )
      stored[[paste(band, endpoint, sep = "|")]] <- list(
        log_ratio = full$log_ratio,
        deleted = deleted
      )
      support[[length(support) + 1L]] <- data.table::rbindlist(
        lapply(triplets, function(x) data.table::data.table(
          band = band,
          endpoint = endpoint,
          year = x$year,
          season = x$season,
          support_size = length(x$support),
          informative = length(x$support) > 1L
        ))
      )
    }
    pedestrian <- stored[[paste(band, PEDESTRIAN, sep = "|")]]
    other <- stored[[paste(band, OTHER_FATAL, sep = "|")]]
    direct <- jackknife_interval(
      pedestrian$log_ratio - other$log_ratio,
      pedestrian$deleted - other$deleted
    )
    estimates[[length(estimates) + 1L]] <- data.table::data.table(
      band = band,
      endpoint = "direct",
      ratio = direct$ratio,
      log_ratio = direct$log_ratio,
      se = direct$se,
      low = direct$low,
      high = direct$high,
      p_value = direct$p_value,
      df = direct$df,
      informative_triplets = NA_integer_,
      total_triplets = NA_integer_,
      degenerate_triplets = NA_integer_
    )
  }
  list(
    estimates = data.table::rbindlist(estimates, fill = TRUE),
    support = data.table::rbindlist(support, fill = TRUE)
  )
}

common_time_analysis <- function(counts) {
  years <- 2001:2024
  estimates <- list()
  decompositions <- list()
  supports <- list()
  for (band in BANDS) {
    prepared <- standardized_inputs(counts, band)
    full <- standardized_estimate(prepared$frame, prepared$years)
    deleted <- lapply(
      prepared$years,
      function(year) standardized_estimate(
        prepared$frame,
        setdiff(prepared$years, year)
      )$estimate
    )
    for (endpoint in c(ENDPOINTS, "direct")) {
      interval <- jackknife_interval(
        full$estimate[[endpoint]],
        vapply(deleted, function(x) x[[endpoint]], numeric(1))
      )
      estimates[[length(estimates) + 1L]] <- data.table::data.table(
        band = band,
        endpoint = endpoint,
        ratio = interval$ratio,
        log_ratio = interval$log_ratio,
        se = interval$se,
        low = interval$low,
        high = interval$high,
        p_value = interval$p_value,
        df = interval$df,
        year_seasons = 2L * length(prepared$years),
        excluded_years = paste(prepared$excluded, collapse = ";")
      )
    }
    decomposition <- data.table::copy(full$parts)
    decomposition[, band := band]
    decompositions[[length(decompositions) + 1L]] <- decomposition
    support <- prepared$frame[, .(
      year,
      season,
      endpoint,
      sham,
      period,
      total,
      fraction,
      gain_side
    )]
    support[, `:=`(band = band, weight = 1 / (2 * length(prepared$years)))]
    supports[[length(supports) + 1L]] <- support
  }
  list(
    estimates = data.table::rbindlist(estimates),
    decomposition = data.table::rbindlist(decompositions),
    support = data.table::rbindlist(supports)
  )
}
