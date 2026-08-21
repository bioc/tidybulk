library(SummarizedExperiment)
library(airway)

test_that("identify_abundant_per_category adds .abundant and respects thresholds", {
  data(airway)
  se <- airway

  expect_warning(
    se_abundant_10 <- tidybulk::identify_abundant_per_category(se, minimum_counts = 10),
    "Neither formula_design and design provided"
  )
  expect_true(".abundant" %in% colnames(rowData(se_abundant_10)))
  expect_type(rowData(se_abundant_10)$.abundant, "logical")
  expect_length(rowData(se_abundant_10)$.abundant, nrow(se))

  se_abundant_100 <- suppressWarnings(
    tidybulk::identify_abundant_per_category(se, minimum_counts = 100)
  )
  expect_true(
    sum(rowData(se_abundant_100)$.abundant) <= sum(rowData(se_abundant_10)$.abundant)
  )
})

test_that("identify_abundant_per_category works with formula_design and minimum_count_per_million", {
  data(airway)
  se <- airway

  se_cpm10 <- tidybulk::identify_abundant_per_category(
    se,
    formula_design = ~ dex,
    minimum_count_per_million = 10
  )
  se_cpm100 <- tidybulk::identify_abundant_per_category(
    se,
    formula_design = ~ dex,
    minimum_count_per_million = 100
  )
  expect_true(
    sum(rowData(se_cpm100)$.abundant) <= sum(rowData(se_cpm10)$.abundant)
  )

  # minimum_count_per_million overrides minimum_counts
  expect_message(
    se_both <- tidybulk::identify_abundant_per_category(
      se,
      formula_design = ~ dex,
      minimum_counts = 100,
      minimum_count_per_million = 100
    ),
    "minimum_count_per_million"
  )
  expect_equal(
    sum(rowData(se_both)$.abundant),
    sum(rowData(se_cpm100)$.abundant)
  )
})

test_that("identify_abundant_per_category works with design and minimum_category", {
  data(airway)
  se <- airway
  design <- model.matrix(~ 0 + dex, data = colData(se))

  se_cat1 <- tidybulk::identify_abundant_per_category(
    se,
    design = design,
    minimum_counts = 10,
    minimum_category = 1
  )
  se_cat2 <- tidybulk::identify_abundant_per_category(
    se,
    design = design,
    minimum_counts = 10,
    minimum_category = 2
  )
  expect_true(
    sum(rowData(se_cat2)$.abundant) <= sum(rowData(se_cat1)$.abundant)
  )
})

test_that("keep_abundant_per_category filters to abundant features", {
  data(airway)
  se <- airway

  se_keep_10 <- tidybulk::keep_abundant_per_category(
    se,
    formula_design = ~ dex,
    minimum_counts = 10
  )
  se_keep_100 <- tidybulk::keep_abundant_per_category(
    se,
    formula_design = ~ dex,
    minimum_counts = 100
  )
  expect_true(nrow(se_keep_100) <= nrow(se_keep_10))
  expect_true(nrow(se_keep_10) <= nrow(se))
  expect_true(all(rowData(se_keep_10)$.abundant))
})

test_that("identify_abundant_per_category validates parameters and design", {
  data(airway)
  se <- airway

  expect_error(
    tidybulk::identify_abundant_per_category(se, formula_design = ~ dex, minimum_counts = -1),
    "minimum_counts"
  )
  expect_error(
    tidybulk::identify_abundant_per_category(se, formula_design = ~ dex, minimum_proportion = 1.5),
    "minimum_proportion"
  )

  # Non-categorical design should error unless coerce_design=TRUE
  design_cont <- model.matrix(~ as.numeric(dex), data = colData(se))
  expect_error(
    tidybulk::identify_abundant_per_category(se, design = design_cont, minimum_counts = 10),
    "design matrix"
  )
  expect_message(
    tidybulk::identify_abundant_per_category(
      se,
      design = design_cont,
      minimum_counts = 10,
      coerce_design = TRUE
    ),
    "coerce"
  )
})

test_that("identify_abundant_per_category force replaces existing .abundant", {
  data(airway)
  se <- airway
  se <- tidybulk::identify_abundant_per_category(se, formula_design = ~ dex, minimum_counts = 10)

  expect_message(
    se2 <- tidybulk::identify_abundant_per_category(
      se,
      formula_design = ~ dex,
      minimum_counts = 100
    ),
    "already exists"
  )
  expect_equal(sum(rowData(se2)$.abundant), sum(rowData(se)$.abundant))

  expect_message(
    se3 <- tidybulk::identify_abundant_per_category(
      se,
      formula_design = ~ dex,
      minimum_counts = 100,
      force = TRUE
    ),
    "force=TRUE"
  )
  expect_true(sum(rowData(se3)$.abundant) <= sum(rowData(se)$.abundant))
})
