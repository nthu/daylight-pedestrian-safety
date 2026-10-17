t_interval <- function(log_ratio, standard_error, degrees_freedom) {
  critical <- stats::qt(0.975, degrees_freedom)
  statistic <- log_ratio / standard_error
  list(
    log_ratio = unname(log_ratio),
    ratio = unname(exp(log_ratio)),
    se = unname(standard_error),
    low = unname(exp(log_ratio - critical * standard_error)),
    high = unname(exp(log_ratio + critical * standard_error)),
    p_value = unname(2 * stats::pt(-abs(statistic), degrees_freedom)),
    df = as.integer(degrees_freedom)
  )
}

jackknife_interval <- function(log_ratio, deleted) {
  deleted <- as.numeric(deleted)
  blocks <- length(deleted)
  if (blocks < 3L || any(!is.finite(deleted))) {
    stop("At least three finite deletion estimates are required.", call. = FALSE)
  }
  variance <- (blocks - 1) / blocks * sum((deleted - mean(deleted))^2)
  standard_error <- sqrt(variance)
  if (!is.finite(standard_error) || standard_error <= 0) {
    stop("Deletion variance is not positive.", call. = FALSE)
  }
  c(t_interval(log_ratio, standard_error, blocks - 1L), list(blocks = blocks))
}

holm_adjust <- function(p_values) stats::p.adjust(p_values, method = "holm")
