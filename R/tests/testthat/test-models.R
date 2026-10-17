testthat::test_that("primary allocation estimates match numerical references", {
  counts <- national_counts()
  reference <- as.data.frame(read_table("national_models.csv", REFERENCE_DIR))
  for (band_value in BANDS) {
    actual <- estimate_allocation(counts, band = band_value)
    expected <- reference[
      reference$label == paste0("primary_", band_value),
      ,
      drop = FALSE
    ]
    expect_fields_close(actual, expected, c("ratio", "low", "high", "se", "df"))
  }
})

testthat::test_that("direct road-user contrasts match numerical references", {
  counts <- national_counts()
  reference <- as.data.frame(read_table("direct_contrasts.csv", REFERENCE_DIR))
  for (band_value in BANDS) {
    pedestrian <- estimate_allocation(counts, band = band_value, endpoint = PEDESTRIAN)
    other <- estimate_allocation(counts, band = band_value, endpoint = OTHER_FATAL)
    actual <- contrast_allocation(
      pedestrian,
      other,
      paste0("pedestrian_vs_other_", band_value)
    )
    expected <- reference[
      reference$contrast == paste0("pedestrian_vs_other_", band_value),
      ,
      drop = FALSE
    ]
    expect_fields_close(actual, expected, c("ratio", "low", "high", "se", "df"))
  }
})

testthat::test_that("fixed-margin conditional estimates match references", {
  tables <- primary_tables()
  reference <- as.data.frame(read_table("conditional_blocks.csv", REFERENCE_DIR))
  for (band_value in BANDS) {
    actual <- exact_block_intervals(tables, band_value, block_sizes = 1)
    expected <- reference[
      reference$band == band_value & reference$block_size_years == 1,
      ,
      drop = FALSE
    ]
    expect_fields_close(actual, expected, c("ratio", "low", "high", "se", "df"))
  }
})

testthat::test_that("year-specific-slope conditional estimates match references", {
  actual <- as.data.frame(triplet_analysis(national_counts())$estimates)
  reference <- as.data.frame(read_table("annual_slope_conditional.csv", REFERENCE_DIR))
  for (band_value in BANDS) {
    for (endpoint_value in c(ENDPOINTS, "direct")) {
      observed <- actual[
        actual$band == band_value & actual$endpoint == endpoint_value,
        ,
        drop = FALSE
      ]
      expected <- reference[
        reference$band == band_value & reference$endpoint == endpoint_value,
        ,
        drop = FALSE
      ]
      expect_fields_close(observed, expected, c("ratio", "low", "high", "se", "df"))
    }
  }
})

testthat::test_that("common-time estimates match references", {
  actual <- as.data.frame(common_time_analysis(national_counts())$estimates)
  reference <- as.data.frame(read_table("common_time_standardized.csv", REFERENCE_DIR))
  for (band_value in BANDS) {
    for (endpoint_value in c(ENDPOINTS, "direct")) {
      observed <- actual[
        actual$band == band_value & actual$endpoint == endpoint_value,
        ,
        drop = FALSE
      ]
      expected <- reference[
        reference$band == band_value & reference$endpoint == endpoint_value,
        ,
        drop = FALSE
      ]
      expect_fields_close(observed, expected, c("ratio", "low", "high", "se", "df"))
    }
  }
})

testthat::test_that("fixed-state and Mississippi estimates match references", {
  state <- as.data.frame(state_analysis(state_counts())$state_effects)
  state_reference <- as.data.frame(read_table("state_stratified.csv", REFERENCE_DIR))
  for (band_value in BANDS) {
    observed <- state[
      state$band == band_value & state$endpoint == PEDESTRIAN,
      ,
      drop = FALSE
    ]
    expected <- state_reference[
      state_reference$band == band_value & state_reference$endpoint == PEDESTRIAN,
      ,
      drop = FALSE
    ]
    expect_fields_close(observed, expected, c("ratio", "low", "high"), tolerance = 2e-6)
  }

  mississippi <- as.data.frame(mississippi_analysis(mississippi_counts())$mississippi_effects)
  reference <- as.data.frame(read_table("mississippi_stratified.csv", REFERENCE_DIR))
  for (band_value in BANDS) {
    observed <- mississippi[
      mississippi$band == band_value & mississippi$endpoint == "all_crashes",
      ,
      drop = FALSE
    ]
    expected <- reference[
      reference$band == band_value & reference$endpoint == "all_crashes",
      ,
      drop = FALSE
    ]
    expect_fields_close(observed, expected, c("ratio", "low", "high"), tolerance = 2e-6)
  }
})

testthat::test_that("legacy event-week coefficients match references", {
  actual <- as.data.frame(estimate_legacy_weeks()$estimates)
  reference <- as.data.frame(read_table("legacy_weekly_estimates.csv", REFERENCE_DIR))
  for (band_value in BANDS) {
    terms <- unique(reference$term[reference$band == band_value])
    for (term_value in terms) {
      observed <- actual[
        actual$band == band_value & actual$term == term_value,
        ,
        drop = FALSE
      ]
      expected <- reference[
        reference$band == band_value & reference$term == term_value,
        ,
        drop = FALSE
      ]
      expect_fields_close(observed, expected, c("ratio", "low", "high", "se", "df"))
    }
  }
})

testthat::test_that("five-center curvature models match references", {
  actual <- as.data.frame(five_center_analysis(full = FALSE)$five_center_curvature)
  reference <- as.data.frame(read_table("five_center_curvature.csv", REFERENCE_DIR))
  reference <- reference[reference$model == "five_center_quadratic", , drop = FALSE]
  for (band_value in BANDS) {
    for (endpoint_value in c(ENDPOINTS, "direct")) {
      observed <- actual[
        actual$model == "five_center_quadratic" &
          actual$band == band_value &
          actual$endpoint == endpoint_value,
        ,
        drop = FALSE
      ]
      expected <- reference[
        reference$band == band_value & reference$endpoint == endpoint_value,
        ,
        drop = FALSE
      ]
      expect_fields_close(observed, expected, c("ratio", "low", "high", "se", "df"))
    }
  }
})

testthat::test_that("exact policy-year multiplier calibration matches references", {
  actual <- as.data.frame(wild_cluster_calibration())
  reference <- as.data.frame(read_table("wild_cluster_calibration.csv", REFERENCE_DIR))
  fields <- c(
    "ratio",
    "wild_p_value",
    "wild_critical_95",
    "wild_low",
    "wild_high",
    "clusters",
    "unique_sign_patterns"
  )
  for (comparison_value in reference$comparison) {
    observed <- actual[actual$comparison == comparison_value, , drop = FALSE]
    expected <- reference[reference$comparison == comparison_value, , drop = FALSE]
    expect_fields_close(observed, expected, fields, tolerance = 1e-6)
  }
})

testthat::test_that("season and policy-era diagnostics match references", {
  season <- as.data.frame(seasonal_heterogeneity())
  season_reference <- as.data.frame(read_table("season_heterogeneity.csv", REFERENCE_DIR))
  for (index in seq_len(nrow(season_reference))) {
    expected <- season_reference[index, , drop = FALSE]
    observed <- season[
      season$band == expected$band &
        season$result_type == expected$result_type &
        season$label == expected$label,
      ,
      drop = FALSE
    ]
    expect_fields_close(observed, expected, c("ratio", "low", "high", "se", "p_value", "df"))
  }

  era <- as.data.frame(policy_era_interaction())
  era_reference <- as.data.frame(read_table("policy_era_interaction.csv", REFERENCE_DIR))
  for (index in seq_len(nrow(era_reference))) {
    expected <- era_reference[index, , drop = FALSE]
    observed <- era[
      era$band == expected$band &
        era$result_type == expected$result_type &
        era$label == expected$label,
      ,
      drop = FALSE
    ]
    expect_fields_close(observed, expected, c("ratio", "low", "high", "se", "p_value", "df"))
  }
})


testthat::test_that("Alaska inclusion sensitivity matches references", {
  actual <- as.data.frame(alaska_inclusion_sensitivity())
  reference <- as.data.frame(read_table("alaska_inclusion.csv", REFERENCE_DIR))
  for (index in seq_len(nrow(reference))) {
    expected <- reference[index, , drop = FALSE]
    observed <- actual[
      actual$geography == expected$geography &
        actual$endpoint == expected$endpoint &
        actual$band == expected$band,
      ,
      drop = FALSE
    ]
    expect_fields_close(observed, expected, c("ratio", "low", "high", "se", "df"))
  }
})
