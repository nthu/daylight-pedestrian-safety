required <- c(
  "data.table",
  "ggplot2",
  "jsonlite",
  "digest",
  "testthat",
  "patchwork",
  "scales"
)

installed <- rownames(installed.packages())
missing <- setdiff(required, installed)
if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  message("All required R packages are already installed.")
}
