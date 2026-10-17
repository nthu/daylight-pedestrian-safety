source(file.path("R", "00_setup.R"))

run_main <- function() {
  counts <- national_counts()
  paths <- result_paths()
  fitted <- list()
  rows <- list()

  for (specification in national_specifications()) {
    name <- specification$name
    arguments <- c(list(counts = counts), specification[names(specification) != "name"])
    fit <- do.call(estimate_allocation, arguments)
    fitted[[name]] <- fit
    rows[[length(rows) + 1L]] <- allocation_summary_row(
      name,
      fit,
      list(
        band = specification$band,
        reference = specification$reference %||% "night",
        endpoint = specification$endpoint %||% PEDESTRIAN,
        geography = specification$geography %||% "DST_main",
        width_minutes = specification$width_minutes %||% 40,
        corrected = specification$corrected %||% TRUE
      )
    )
  }

  comparisons <- list(
    c("primary_evening", "primary_morning", "evening_vs_morning"),
    c("primary_evening", "other_fatal_crashes_evening", "pedestrian_vs_other_evening"),
    c("primary_morning", "other_fatal_crashes_morning", "pedestrian_vs_other_morning"),
    c("early_evening", "late_evening", "early_vs_late_evening"),
    c("east_evening", "west_evening", "east_vs_west_evening"),
    c("early_morning", "late_morning", "early_vs_late_morning"),
    c("east_morning", "west_morning", "east_vs_west_morning")
  )
  direct <- data.table::rbindlist(lapply(comparisons, function(item) {
    data.table::as.data.table(contrast_allocation(
      fitted[[item[1]]],
      fitted[[item[2]]],
      item[3]
    ))
  }), fill = TRUE)
  estimates <- data.table::rbindlist(rows, fill = TRUE)

  write_table(estimates, file.path(paths$estimates, "national_models.csv"))
  write_table(direct, file.path(paths$estimates, "direct_contrasts.csv"))

  cat("\nNational primary estimates\n")
  print(estimates[specification %chin% c("primary_evening", "primary_morning"),
    .(specification, ratio, low, high)])
  cat("\nDirect road-user comparisons\n")
  print(direct[, .(contrast, ratio, low, high)])
  invisible(list(models = fitted, estimates = estimates, direct = direct))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

if (sys.nframe() == 0L) run_main()
