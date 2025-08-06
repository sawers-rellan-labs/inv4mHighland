test_that("Package structure is correct", {
  expect_true(file.exists("R/clayton_spatial_analysis.R"))
  expect_true(file.exists("R/spatial_models.R"))
  expect_true(file.exists("R/variogram_analysis.R"))
  expect_true(file.exists("DESCRIPTION"))
  expect_true(file.exists("NAMESPACE"))
})

test_that("Main functions exist", {
  expect_true(exists("run_clayton_spatial_analysis"))
  expect_true(exists("fit_all_models"))
  expect_true(exists("calculate_scaled_variogram"))
  expect_true(exists("create_spatial_plot"))
})

test_that("Utility functions work", {
  # Test safe_fixef with a simple model
  test_data <- data.frame(
    y = rnorm(20),
    x = rnorm(20)
  )
  
  test_model <- lm(y ~ x, data = test_data)
  effects <- safe_fixef(test_model)
  
  expect_true(!is.null(effects))
  expect_true(is.numeric(effects))
  expect_equal(length(effects), 2)  # intercept + slope
})