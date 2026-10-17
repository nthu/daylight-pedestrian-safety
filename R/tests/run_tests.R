repository <- normalizePath(getwd(), mustWork = TRUE)
if (!dir.exists(file.path(repository, "R", "tests", "testthat"))) {
  stop("Run this script from the repository root.", call. = FALSE)
}

Sys.setenv(DAYLIGHT_REPO_ROOT = repository)

result <- testthat::test_dir(
  file.path(repository, "R", "tests", "testthat"),
  reporter = "summary",
  stop_on_failure = TRUE
)

invisible(result)
