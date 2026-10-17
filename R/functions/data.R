find_repository_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(current, "data")) &&
        dir.exists(file.path(current, "R")) &&
        dir.exists(file.path(current, "python"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the repository root.", call. = FALSE)
    }
    current <- parent
  }
}

REPOSITORY <- find_repository_root()
DATA_DIR <- file.path(REPOSITORY, "data")
RESULTS_DIR <- file.path(REPOSITORY, "results", "R")

read_table <- function(name, folder = DATA_DIR) {
  path <- file.path(folder, name)
  if (!file.exists(path)) {
    stop("Missing input: ", path, call. = FALSE)
  }
  data.table::fread(path, showProgress = FALSE)
}

national_counts <- function() read_table("national.csv.gz")
state_counts <- function() read_table("states.csv.gz")
mississippi_counts <- function() read_table("mississippi.csv.gz")
reference_counts <- function() read_table("reference_holiday.csv.gz")
traffic_counts <- function() read_table("traffic.csv")
primary_tables <- function() read_table("primary_tables.csv")
policy_calendar <- function() read_table("calendar.csv")
sensitivity_table <- function(name) read_table(name, file.path(DATA_DIR, "sensitivity"))
support_table <- function(name) read_table(name, file.path(DATA_DIR, "support"))

result_paths <- function() {
  paths <- list(
    root = RESULTS_DIR,
    estimates = file.path(RESULTS_DIR, "estimates"),
    tables = file.path(RESULTS_DIR, "tables"),
    figures = file.path(RESULTS_DIR, "figures")
  )
  invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))
  paths
}

write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(data.table::as.data.table(x), path, na = "", dateTimeAs = "write.csv")
  invisible(path)
}
