arguments <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% arguments

for (script in sprintf("R/%02d_%s.R", 1:5, c(
  "main",
  "sensitivity",
  "registry",
  "traffic",
  "outputs"
))) {
  source(script)
}

stages <- list(
  "01_main.R" = function() run_main(),
  "02_sensitivity.R" = function() run_sensitivity(full = !quick),
  "03_registry.R" = function() run_registry(),
  "04_traffic.R" = function() run_traffic(),
  "05_outputs.R" = function() run_outputs()
)

for (name in names(stages)) {
  cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
  cat(name, "\n")
  cat(paste(rep("=", 72), collapse = ""), "\n", sep = "")
  stages[[name]]()
}
