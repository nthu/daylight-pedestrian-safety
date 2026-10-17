source(file.path("R", "00_setup.R"))

run_registry <- function() {
  paths <- result_paths()
  state <- state_analysis(state_counts())
  mississippi <- mississippi_analysis(mississippi_counts())
  outputs <- c(state, mississippi)
  for (name in names(outputs)) {
    write_table(outputs[[name]], file.path(paths$estimates, paste0(name, ".csv")))
  }

  cat("\nFixed-state pedestrian estimates\n")
  print(state$state_effects[endpoint == PEDESTRIAN, .(band, ratio, low, high)])
  cat("\nMississippi all-crash estimates\n")
  print(mississippi$mississippi_effects[
    endpoint == "all_crashes",
    .(band, ratio, low, high, holm_p)
  ])
  invisible(outputs)
}

if (sys.nframe() == 0L) run_registry()
