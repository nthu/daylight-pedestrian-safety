repository <- Sys.getenv("DAYLIGHT_REPO_ROOT")
if (!nzchar(repository)) stop("DAYLIGHT_REPO_ROOT is not set.", call. = FALSE)
source(file.path(repository, "R", "00_setup.R"))

REFERENCE_DIR <- file.path(REPOSITORY, "R", "reference")

expect_fields_close <- function(actual, expected, fields, tolerance = 3e-6) {
  for (field in fields) {
    testthat::expect_equal(
      as.numeric(actual[[field]]),
      as.numeric(expected[[field]]),
      tolerance = tolerance,
      info = field
    )
  }
}
