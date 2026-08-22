#' Add edgeR dispersion estimates to rowData
#'
#' Port of HPCell's `se_estimate_dispersion()`
#' (<https://github.com/MangiolaLaboratory/HPCell/blob/ed763c6f53cfedf63d3a3353ec5ffcfd09bfe0e6/R/differential_expression.R#L36>):
#' tagwise dispersion when there are fewer than 1000 samples, otherwise
#' trended dispersion on a random subset of up to 2000 samples. The design is
#' a fixed-effect formula: edgeR has no random effects, so write the analogue
#' yourself (`~ dex + cell` for `~ dex + (1 | cell)`, `~ dex * cell` for
#' `~ dex + (dex | cell)`). [estimate_gene()] does not compute this;
#' [estimate()] calls this helper once on the full object before iterating
#' genes.
#'
#' @param .data A `SummarizedExperiment`
#' @param formula_abundance Fixed-effect model for the mean, used to build
#'   the edgeR design. This is not the mixed model you pass to
#'   [estimate_gene()]: random-effect terms (`|`) are rejected. Translate
#'   them yourself, e.g. `(1 | donor)` as `donor` and
#'   `(1 + treatment | donor)` as `donor + treatment:donor`.
#' @param abundance Assay name (default `"counts"`).
#' @param dispersion_column Name of the `rowData` column for tagwise (or
#'   trended) dispersion \eqn{\phi_g} (default `"dispersion"`).
#' @param dispersion_degrees_freedom_column Name of the `rowData` column for
#'   the effective degrees of freedom \eqn{d_{eff}} behind that estimate
#'   (default `"dispersion_degrees_freedom"`).
#'
#' @details
#' Two columns are written, named by `dispersion_column` and
#' `dispersion_degrees_freedom_column`. The first holds \eqn{\phi_g}; the
#' second holds the effective degrees of freedom behind that estimate.
#' [estimate_gene()] turns the pair into a **prior** on the negative binomial
#' shape, not a plug-in: \eqn{\phi_g} locates the prior and \eqn{d_{eff}} sets
#' its tightness, so the gene-wise likelihood can still pull the posterior
#' away from the edgeR value when the counts disagree.
#'
#' The effective degrees of freedom are
#' \deqn{d_{eff} = (n - \mathrm{ncol}(design)) + \mathrm{prior.df}}
#' The first term is the residual degrees of freedom of the fixed-effect
#' design; `estimateDisp()` does not return it, so it is computed here. The
#' second is the prior degrees of freedom that `estimateDisp()` does return,
#' quantifying how far each gene-wise estimate is shrunk toward the
#' mean-dispersion trend by empirical Bayes. `estimateTrendedDisp()` (the
#' branch used above 1000 samples) performs no such shrinkage and returns no
#' `prior.df`, so there \eqn{d_{eff}} is the residual degrees of freedom alone.
#'
#' A grouping factor written as a fixed effect with close to one level per
#' sample will exhaust the design and leave no residual degrees of freedom.
#' That is an error rather than a column of `NA`s: pass a simpler
#' `formula_abundance` for this step if it happens.
#'
#' The link from degrees of freedom to a standard deviation runs through the
#' chi-square distribution of the estimator. For \eqn{s^2 \sim \sigma^2
#' \chi^2_d / d}, the log of the estimate has variance
#' \deqn{\mathrm{Var}(\log s^2) = \psi'(d/2)}
#' with \eqn{\psi'} the trigamma function, so
#' \eqn{\mathrm{SD}(\log \hat\phi_g) \approx \sqrt{\psi'(d_{eff}/2)}}. The
#' familiar \eqn{\sqrt{2/d}} is the large-\eqn{d} approximation to this and is
#' 12% too small at \eqn{d = 4}, so the trigamma form is used directly. This
#' is the same empirical-Bayes variance model that limma and edgeR use to
#' moderate gene-wise variances; it is exact for a scaled chi-square and is
#' applied to the negative binomial dispersion by analogy.
#'
#' The same \eqn{d_{eff}} supports a second, conjugate parameterisation, which
#' [estimate_gene()] offers as `shape_prior = "gamma"`. Inverting the scaled
#' inverse chi-square gives a gamma, so the precision \eqn{1/\phi_g} -- the
#' quantity brms calls `shape` -- has prior
#' \deqn{\mathrm{Gamma}(d_{eff}/2,\ \mathrm{rate} = d_{eff}\phi_g/2)}
#' with mean \eqn{1/\phi_g}. Because \eqn{\mathrm{Var}(\log X) = \psi'(a)} for
#' a gamma of shape \eqn{a}, this reproduces \eqn{\psi'(d_{eff}/2)} exactly:
#' the two routes are the same calculation, since \eqn{\chi^2_d} is itself
#' \eqn{\mathrm{Gamma}(d/2, \mathrm{scale} = 2)}. Its coefficient of variation
#' is \eqn{\sqrt{2/d_{eff}}}, which is where that familiar approximation comes
#' from.
#'
#' The two forms are not reparameterisations of one another, and differ in two
#' ways worth knowing.
#'
#' They centre different summaries of the shape on edgeR's estimate. brms
#' gives the shape submodel a log link, so a Student-t on the intercept is
#' symmetric in \eqn{\log(\mathrm{shape})} and places the *median* of the
#' shape at \eqn{1/\phi_g}. The gamma places its *mean* there, and a gamma's
#' median lies below its mean. On the log scale the centres differ by
#' \deqn{\psi(d_{eff}/2) - \log(d_{eff}/2) \approx -1/d_{eff}}
#' with \eqn{\psi} the digamma function: about -0.11, or 0.22 prior standard
#' deviations, at \eqn{d_{eff} = 9.8}, and vanishing as samples accumulate.
#' Neither is wrong. Mean-centring is what the conjugate hierarchy dictates,
#' since a scaled inverse chi-square on \eqn{\phi_g} implies
#' \eqn{E[1/\phi_g] = 1/\phi_g}; median-centring is what an additive offset on
#' a log link implies.
#'
#' They also differ in tail weight, which is why the Student-t is the default.
#' Exponentiating a Student-t leaves a prior on the shape with polynomial
#' tails and no finite mean, proper but very permissive, so a gene whose true
#' dispersion is far from the shrunken edgeR estimate can still escape. The
#' gamma decays exponentially and holds such a gene closer to the trend. The
#' gamma is also left-skewed in \eqn{\log(\mathrm{shape})} (skewness
#' \eqn{\psi''(a) / \psi'(a)^{3/2}}, about -0.47 here) where the Student-t is
#' symmetric.
#'
#' Note that this standard deviation describes how precisely edgeR estimated
#' its own \eqn{\phi_g}. It does not account for the fact that the design here
#' fits grouping factors as fixed rather than partially pooled, nor for a
#' zero-inflated likelihood in which the `zi` component absorbs part of the
#' overdispersion. Treating them as fixed is the conservative choice: it
#' spends the full degrees of freedom that shrinkage would have given back,
#' so \eqn{d_{eff}} understates rather than overstates the information behind
#' \eqn{\phi_g}.
#'
#' @references
#' Smyth GK (2004). Linear models and empirical Bayes methods for assessing
#' differential expression in microarray experiments. *Statistical
#' Applications in Genetics and Molecular Biology* 3(1).
#' [PDF](https://gksmyth.github.io/pubs/ebayes.pdf) — derives the
#' \eqn{\mathrm{Var}(\log s^2) = \psi'(d/2)} result and the prior degrees of
#' freedom used to moderate it.
#'
#' McCarthy DJ, Chen Y, Smyth GK (2012). Differential expression analysis of
#' multifactor RNA-Seq experiments with respect to biological variation.
#' *Nucleic Acids Research* 40(10):4288-4297.
#' [PMC3378882](https://pmc.ncbi.nlm.nih.gov/articles/PMC3378882/) —
#' Cox-Reid adjusted profile likelihood dispersion conditional on a design.
#'
#' Robinson MD, McCarthy DJ, Smyth GK (2010). edgeR: a Bioconductor package
#' for differential expression analysis of digital gene expression data.
#' *Bioinformatics* 26(1):139-140.
#' [PMC2796818](https://pmc.ncbi.nlm.nih.gov/articles/PMC2796818/)
#'
#' Phipson B, Lee S, Majewski IJ, Alexander WS, Smyth GK (2016). Robust
#' hyperparameter estimation protects against hypervariable genes and improves
#' power to detect differential expression. *Annals of Applied Statistics*
#' 10(2):946-963. \doi{10.1214/16-AOAS920} — `robust = TRUE`, which makes
#' `prior.df` gene-specific and possibly infinite.
#'
#' [`estimateDisp()` reference manual](https://rdrr.io/bioc/edgeR/man/estimateDisp.html),
#' [edgeR on Bioconductor](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
#'
#' @return `.data` with `rowData(.data)[[dispersion_column]]` and
#'   `rowData(.data)[[dispersion_degrees_freedom_column]]` filled. The raw
#'   edgeR object is in `metadata(.)$tidybulk$estimateDisp`.
#'
#' @examples
#' \dontrun{
#' data("airway", package = "airway")
#' se <- estimate_dispersion(airway[1:150, ], ~ dex + cell)
#' }
#'
#' @importFrom SummarizedExperiment assay rowData colData "rowData<-"
#' @name estimate_dispersion
#' @docType methods
#' @rdname estimate_dispersion-methods
#' @export
setGeneric("estimate_dispersion", function(.data,
                                           formula_abundance,
                                           abundance = "counts",
                                           dispersion_column = "dispersion",
                                           dispersion_degrees_freedom_column = "dispersion_degrees_freedom")
  standardGeneric("estimate_dispersion"))

.estimate_dispersion_se <- function(.data,
                                    formula_abundance,
                                    abundance = "counts",
                                    dispersion_column = "dispersion",
                                    dispersion_degrees_freedom_column = "dispersion_degrees_freedom") {
  dispersion_column <- check_dispersion_name(dispersion_column)
  dispersion_degrees_freedom_column <-
    check_degrees_freedom_name(dispersion_degrees_freedom_column)
  if (identical(dispersion_column, dispersion_degrees_freedom_column)) {
    stop(
      "tidybulk says: `dispersion_column` and `dispersion_degrees_freedom_column` must name different columns.",
      call. = FALSE
    )
  }
  check_and_install_packages("edgeR")

  n_sample <- ncol(.data)

  check_formula(formula_abundance, force_fixed_effects = TRUE)

# If the sample size is less than 1000, use the tagwise dispersion
  if (n_sample < 1000L) {
    design <- dispersion_design(.data, formula_abundance)
    check_residual_df(design, n_sample)
    counts <- assay(.data, abundance)
    dispersion_object <- edgeR::estimateDisp(counts, design = design)
    disp <- dispersion_object$tagwise.dispersion
    d_eff <- (n_sample - ncol(design)) + dispersion_object$prior.df
  } 
  
  # If the sample size is greater than 1000, use the trended dispersion
  else {
    sampled <- sample(seq_len(n_sample), size = min(n_sample, 2000L))
    se_sub <- .data[, sampled, drop = FALSE]
    design <- dispersion_design(se_sub, formula_abundance)
    check_residual_df(design, ncol(se_sub))
    counts <- assay(se_sub, abundance)
    dispersion_object <- edgeR::estimateTrendedDisp(
      counts,
      design = design,
      subset = 1000,
      rowsum.filter = 10
    )
    disp <- dispersion_object
    d_eff <- ncol(se_sub) - ncol(design)
  }

  rowData(.data)[[dispersion_column]] <- as.numeric(disp)
  rowData(.data)[[dispersion_degrees_freedom_column]] <-
    rep_len(as.numeric(d_eff), nrow(.data))

  .data <- attach_to_metadata(.data, dispersion_object, "estimateDisp")
  .data <- memorise_methods_used(.data, "edger")
  rlang::inform(
    "tidybulk says: to access the raw results do `metadata(.)$tidybulk$estimateDisp`",
    .frequency_id = "Access estimateDisp results",
    .frequency = "always"
  )
  .data
}

#' estimate_dispersion
#'
#' @docType methods
#' @rdname estimate_dispersion-methods
#'
#' @return A `SummarizedExperiment` object
setMethod("estimate_dispersion",
          "SummarizedExperiment",
          .estimate_dispersion_se)

#' estimate_dispersion
#'
#' @docType methods
#' @rdname estimate_dispersion-methods
#'
#' @return A `SummarizedExperiment` object
setMethod("estimate_dispersion",
          "RangedSummarizedExperiment",
          .estimate_dispersion_se)

# A grouping factor with close to one level per sample can exhaust the
# design. edgeR would return NA for every gene, so stop instead of handing
# back a column of NAs.
check_residual_df <- function(design, n_sample) {
  if (n_sample - ncol(design) > 0L) {
    return(invisible(NULL))
  }
  stop(
    "tidybulk says: The dispersion design has ", ncol(design), " coefficients for ",
    n_sample, " samples, leaving no residual degrees of freedom. ",
    "A grouping factor with nearly one level per sample will exhaust the ",
    "design. Pass a simpler `formula_abundance` to estimate_dispersion().",
    call. = FALSE
  )
}

dispersion_design <- function(se, formula) {
  stats::model.matrix(
    formula,
    data = droplevels(as.data.frame(SummarizedExperiment::colData(se)))
  )
}

check_dispersion_name <- function(x) x
check_degrees_freedom_name <- function(x) x
