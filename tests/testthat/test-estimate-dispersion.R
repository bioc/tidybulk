airway_se <- function(n_genes = 150) {
  skip_if_not_installed("airway")
  data("airway", package = "airway")
  airway[seq_len(n_genes), ]
}

simulated_se <- function(n_genes = 100) {
  set.seed(1)
  n_donor <- 6L
  n_rep <- 2L
  n <- n_donor * 2L * n_rep
  donor <- factor(rep(paste0("d", seq_len(n_donor)), each = 2L * n_rep))
  treatment <- factor(rep(rep(c("ctrl", "trt"), each = n_rep), n_donor))
  true_dispersion <- exp(stats::rnorm(n_genes, mean = log(0.1), sd = 0.4))
  mu <- 80 * (1 + 0.4 * as.integer(treatment == "trt"))
  counts <- matrix(
    stats::rnbinom(n_genes * n, mu = rep(mu, each = n_genes), size = 1 / true_dispersion),
    nrow = n_genes
  )
  colnames(counts) <- paste0("s", seq_len(n))
  rownames(counts) <- paste0("g", seq_len(n_genes))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(treatment = treatment, donor = donor),
    rowData = S4Vectors::DataFrame(true_dispersion = true_dispersion)
  )
}

test_that("estimate_dispersion writes a positive tagwise column on airway", {
  skip_if_not_installed("edgeR")

  se <- airway_se(n_genes = 150)
  out <- suppressMessages(estimate_dispersion(se, ~ dex + cell, abundance = "counts"))
  rd <- SummarizedExperiment::rowData(out)
  expect_equal(length(rd$dispersion_shrinked), nrow(se))
  expect_true(all(is.finite(rd$dispersion_shrinked)))
  expect_true(all(rd$dispersion_shrinked > 0))
  expect_equal(length(rd$dispersion_trended), nrow(se))
  expect_true(all(is.finite(rd$dispersion_trended)))
  expect_true(all(rd$dispersion_trended > 0))
})

test_that("estimate_dispersion records d_eff = df.residual + prior.df", {
  skip_if_not_installed("edgeR")

  se <- airway_se(n_genes = 150)
  out <- suppressMessages(estimate_dispersion(se, ~ dex + cell, abundance = "counts"))
  d_eff <- SummarizedExperiment::rowData(out)$dispersion_degrees_freedom
  expect_equal(length(d_eff), nrow(se))

  design <- stats::model.matrix(
    ~ dex + cell,
    data = droplevels(as.data.frame(SummarizedExperiment::colData(se)))
  )
  fit <- edgeR::estimateDisp(
    SummarizedExperiment::assay(se, "counts"),
    design = design
  )
  expect_equal(
    unique(d_eff),
    (ncol(se) - ncol(design)) + fit$prior.df
  )
})

test_that("estimate_dispersion stores the edgeR object in tidybulk metadata", {
  skip_if_not_installed("edgeR")

  se <- airway_se(n_genes = 150)
  expect_message(
    out <- estimate_dispersion(se, ~ dex + cell, abundance = "counts"),
    "metadata\\(\\.\\)\\$tidybulk\\$estimateDisp"
  )
  tb <- S4Vectors::metadata(out)$tidybulk
  expect_true("estimateDisp" %in% names(tb))
  expect_equal(
    as.numeric(tb$estimateDisp$tagwise.dispersion),
    as.numeric(SummarizedExperiment::rowData(out)$dispersion_shrinked)
  )
  expect_equal(
    as.numeric(tb$estimateDisp$trended.dispersion),
    as.numeric(SummarizedExperiment::rowData(out)$dispersion_trended)
  )
  expect_true("edger" %in% tb$methods_used)
})

test_that("estimate_dispersion output column names are arguments", {
  skip_if_not_installed("edgeR")

  se <- airway_se(n_genes = 150)
  out <- suppressMessages(estimate_dispersion(
    se,
    ~ dex + cell,
    abundance = "counts",
    dispersion_column = "phi",
    trended_dispersion_column = "phi_trend",
    dispersion_degrees_freedom_column = "phi_deff"
  ))
  rd <- SummarizedExperiment::rowData(out)
  expect_true(all(c("phi", "phi_trend", "phi_deff") %in% names(rd)))
  expect_false(any(c("dispersion_shrinked", "dispersion_trended", "dispersion_degrees_freedom") %in% names(rd)))
  expect_true(all(rd$phi > 0))
  expect_true(all(rd$phi_trend > 0))
  expect_length(unique(rd$phi_deff), 1L)

  expect_error(
    suppressMessages(estimate_dispersion(
      se, ~dex,
      dispersion_column = "x",
      trended_dispersion_column = "x"
    )),
    "different columns"
  )
})

test_that("estimate_dispersion rejects random-effect terms", {
  skip_if_not_installed("edgeR")
  se <- simulated_se()
  expect_error(
    estimate_dispersion(se, ~ treatment + (1 | donor)),
    "cannot include random effects"
  )
  expect_error(
    estimate_dispersion(se, ~ treatment + (treatment | donor)),
    "cannot include random effects"
  )
})

test_that("estimate_dispersion fails when the design is saturated", {
  skip_if_not_installed("edgeR")

  se <- simulated_se()
  se$sid <- factor(colnames(se))
  expect_error(
    estimate_dispersion(se, ~ treatment + sid),
    "no residual degrees of freedom"
  )
})

test_that("check_formula rejects a non-formula", {
  skip_if_not_installed("edgeR")
  expect_error(
    estimate_dispersion(simulated_se(), "treatment"),
    "must be an R formula"
  )
})

test_that("estimate_dispersion spends degrees of freedom on the design", {
  skip_if_not_installed("edgeR")

  se <- simulated_se()
  out <- suppressMessages(estimate_dispersion(se, ~ treatment * donor))
  disp <- SummarizedExperiment::rowData(out)$dispersion_shrinked
  d_eff <- SummarizedExperiment::rowData(out)$dispersion_degrees_freedom
  expect_true(all(is.finite(disp) & disp > 0))

  design <- stats::model.matrix(
    ~ treatment * donor,
    data = as.data.frame(SummarizedExperiment::colData(se))
  )
  expect_equal(ncol(design), 12L)
  fit <- edgeR::estimateDisp(
    SummarizedExperiment::assay(se, "counts"),
    design = design
  )
  expect_equal(unique(d_eff), (ncol(se) - ncol(design)) + fit$prior.df)

  d_eff_additive <- SummarizedExperiment::rowData(
    suppressMessages(estimate_dispersion(se, ~ treatment + donor))
  )$dispersion_degrees_freedom
  expect_gt(unique(d_eff_additive), unique(d_eff))
})

test_that("estimate_dispersion recovers the dispersion it was simulated with", {
  skip_if_not_installed("edgeR")

  se <- simulated_se()
  out <- suppressMessages(estimate_dispersion(se, ~ treatment + donor))
  expect_gt(
    stats::cor(
      SummarizedExperiment::rowData(out)$dispersion_shrinked,
      SummarizedExperiment::rowData(out)$true_dispersion,
      method = "spearman"
    ),
    0.6
  )
})

test_that("estimate_dispersion is an SE method", {
  expect_error(
    estimate_dispersion(mtcars, ~ dex),
    "unable to find an inherited method"
  )
})
