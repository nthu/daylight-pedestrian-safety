options(stringsAsFactors = FALSE, warn = 1)

setup_repository <- Sys.getenv("DAYLIGHT_REPO_ROOT")
if (!nzchar(setup_repository)) {
  current <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(current, "data")) &&
        dir.exists(file.path(current, "R")) &&
        dir.exists(file.path(current, "python"))) {
      setup_repository <- current
      break
    }
    parent <- dirname(current)
    if (identical(parent, current)) stop("Could not locate the repository root.", call. = FALSE)
    current <- parent
  }
}

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

required_packages <- c(
  "data.table",
  "ggplot2",
  "jsonlite",
  "digest",
  "testthat",
  "patchwork",
  "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the required R packages before running the analysis: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

source(file.path(setup_repository, "R", "functions", "data.R"))
source(file.path(setup_repository, "R", "functions", "inference.R"))
source(file.path(setup_repository, "R", "functions", "allocation.R"))
source(file.path(setup_repository, "R", "functions", "exact.R"))
source(file.path(setup_repository, "R", "functions", "targeted.R"))
source(file.path(setup_repository, "R", "functions", "registry.R"))
source(file.path(setup_repository, "R", "functions", "legacy.R"))
source(file.path(setup_repository, "R", "functions", "specifications.R"))
source(file.path(setup_repository, "R", "functions", "curvature.R"))
source(file.path(setup_repository, "R", "functions", "review_diagnostics.R"))
source(file.path(setup_repository, "R", "functions", "plots.R"))

BANDS <- c("evening", "morning")
SEASONS <- c("spring", "autumn")
PEDESTRIAN <- "fatal_pedestrian_crashes"
OTHER_FATAL <- "other_fatal_crashes"
ENDPOINTS <- c(PEDESTRIAN, OTHER_FATAL)
