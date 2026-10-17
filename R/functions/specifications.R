national_specifications <- function() {
  rows <- list()
  for (band in BANDS) {
    rows <- c(rows, list(
      list(name = paste0("primary_", band), band = band),
      list(name = paste0("fixed_clock_", band), band = paste0(band, "_fixed")),
      list(name = paste0("actual_only_", band), band = band, corrected = FALSE),
      list(name = paste0("all_fatal_crashes_", band), band = band, endpoint = "all_fatal_crashes"),
      list(name = paste0("other_fatal_crashes_", band), band = band, endpoint = OTHER_FATAL),
      list(name = paste0("early_", band), band = band, years = 2001:2010),
      list(name = paste0("late_", band), band = band, years = 2011:2024),
      list(
        name = paste0("without_COVID_", band),
        band = band,
        years = setdiff(2001:2024, c(2020, 2021))
      ),
      list(name = paste0("east_", band), band = band, geography = "DST_east"),
      list(name = paste0("west_", band), band = band, geography = "DST_west"),
      list(name = paste0("inner_", band), band = band, width_minutes = 20),
      list(name = paste0("day_reference_", band), band = band, reference = "noon"),
      list(name = paste0("spring_", band), band = band, seasons = "spring"),
      list(name = paste0("autumn_", band), band = band, seasons = "autumn")
    ))
  }
  rows
}

calendar_specifications <- c(
  matched_14d_gap35 = "calendar_counts_matched_14d_gap35.csv",
  nonoverlap_14d_gap42 = "calendar_counts_nonoverlap_14d_gap42.csv",
  short_7d_gap35 = "calendar_counts_short_7d_gap35.csv",
  long_21d_gap35 = "calendar_counts_long_21d_gap35.csv"
)

allocation_summary_row <- function(name, fit, metadata = list()) {
  data.table::as.data.table(c(
    list(specification = name),
    metadata,
    list(
      ratio = fit$ratio,
      low = fit$low,
      high = fit$high,
      se = fit$se,
      df = fit$df,
      p_value = fit$p_value
    )
  ))
}
