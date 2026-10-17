testthat::test_that("core aggregate tables have the expected support", {
  testthat::expect_equal(nrow(national_counts()), 48384L)
  testthat::expect_equal(nrow(state_counts()), 446688L)
  testthat::expect_equal(nrow(mississippi_counts()), 516600L)
  testthat::expect_equal(nrow(traffic_counts()), 352L)

  national <- national_counts()
  testthat::expect_true(all(national$events >= 0))
  testthat::expect_equal(sort(unique(national$year)), 2001:2024)
  testthat::expect_equal(sort(unique(national$sham)), c(-35, 0, 35))
  testthat::expect_setequal(unique(national$period), c("pre", "post"))

  states <- state_counts()
  testthat::expect_equal(data.table::uniqueN(states$state), 47L)
  testthat::expect_true(all(states$events >= 0))

  mississippi <- mississippi_counts()
  testthat::expect_equal(data.table::uniqueN(mississippi$county_key), 82L)
  testthat::expect_equal(sort(unique(mississippi$year)), 2004:2024)
})

testthat::test_that("policy windows are equal-duration and weekday matched", {
  calendar <- policy_calendar()
  testthat::expect_equal(nrow(calendar), 24L * 2L * 3L * 2L)
  testthat::expect_true(all(calendar$days == 14L))
  calendar[, start_date := as.Date(start)]
  calendar[, end_date := as.Date(end)]
  testthat::expect_true(all(calendar$end_date - calendar$start_date == 13))
  paired <- data.table::dcast(
    calendar,
    year + season + sham ~ period,
    value.var = "start_date"
  )
  testthat::expect_true(all(as.POSIXlt(paired$pre)$wday == as.POSIXlt(paired$post)$wday))
})

testthat::test_that("released aggregate tables exclude event identifiers", {
  forbidden <- c(
    "event_id", "st_case", "veh_no", "per_no", "vin", "ssn",
    "person_id", "local_datetime", "longitude", "station_id"
  )
  paths <- list.files(DATA_DIR, recursive = TRUE, full.names = TRUE)
  paths <- paths[grepl("\\.csv(?:\\.gz)?$", paths)]
  for (path in paths) {
    header <- names(data.table::fread(path, nrows = 0, showProgress = FALSE))
    lower <- tolower(header)
    hits <- intersect(lower, forbidden)
    if (basename(path) == "solar_grid.csv") {
      hits <- setdiff(hits, "latitude")
    }
    testthat::expect_equal(length(hits), 0L, info = path)
  }
})

testthat::test_that("data manifest matches the released files", {
  lines <- readLines(file.path(DATA_DIR, "manifest.sha256"), warn = FALSE)
  testthat::expect_true(length(lines) > 0L)
  for (line in lines) {
    parts <- strsplit(line, "  ", fixed = TRUE)[[1]]
    testthat::expect_equal(length(parts), 2L)
    path <- file.path(DATA_DIR, parts[2])
    testthat::expect_true(file.exists(path), info = parts[2])
    actual <- digest::digest(file = path, algo = "sha256", serialize = FALSE)
    testthat::expect_identical(actual, parts[1])
  }
})
