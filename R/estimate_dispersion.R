# Internal. edgeR tagwise (or trended if n >= 1000) written to rowData.
# `(1 | g)` becomes `g` so grouping spends residual df rather than inflating phi.

fixed_effects_formula <- function(formula) {
  labels <- attr(stats::terms(formula), "term.labels")
  is_re <- grepl("|", labels, fixed = TRUE)
  if (!any(is_re)) {
    return(formula)
  }
  extra <- unlist(lapply(labels[is_re], function(term) {
    sides <- trimws(strsplit(term, "|", fixed = TRUE)[[1]])
    lhs <- stats::terms(stats::as.formula(paste("~", sides[[1]])))
    group <- attr(stats::terms(stats::as.formula(paste("~", sides[[length(sides)]]))), "term.labels")
    out <- if (attr(lhs, "intercept") == 1L) group else character()
    lhs_lab <- attr(lhs, "term.labels")
    if (length(lhs_lab) && length(group)) {
      out <- c(out, as.vector(outer(lhs_lab, group, paste, sep = ":")))
    }
    out
  }), use.names = FALSE)
  stats::reformulate(unique(c(labels[!is_re], extra)))
}

estimate_dispersion <- function(.data, formula_abundance, abundance = "counts") {
  check_and_install_packages("edgeR")
  check_formula(formula_abundance)
  if (grepl("|", paste(deparse(formula_abundance), collapse = ""), fixed = TRUE)) {
    converted <- fixed_effects_formula(formula_abundance)
    message(
      "tidybulk says: ", paste(deparse(formula_abundance), collapse = ""),
      " converted to ", paste(deparse(converted), collapse = ""),
      " for edgeR dispersion (edgeR has no random effects)."
    )
    formula_abundance <- converted
  }

  n_sample <- ncol(.data)
  design <- stats::model.matrix(
    formula_abundance,
    data = droplevels(as.data.frame(colData(.data)))
  )
  if (n_sample <= ncol(design)) {
    stop(
      "tidybulk says: The dispersion design has ", ncol(design),
      " coefficients for ", n_sample,
      " samples, leaving no residual degrees of freedom.",
      call. = FALSE
    )
  }

  if (n_sample < 1000L) {
    counts <- assay(.data, abundance)
    fit <- edgeR::estimateDisp(counts, design = design)
    disp <- fit$tagwise.dispersion
    trended <- fit$trended.dispersion
    d_eff <- (n_sample - ncol(design)) + fit$prior.df
  } else {
    sampled <- sample(seq_len(n_sample), size = min(n_sample, 2000L))
    se_sub <- .data[, sampled, drop = FALSE]
    design <- stats::model.matrix(
      formula_abundance,
      data = droplevels(as.data.frame(colData(se_sub)))
    )
    fit <- edgeR::estimateTrendedDisp(
      assay(se_sub, abundance),
      design = design,
      subset = 1000,
      rowsum.filter = 10
    )
    disp <- trended <- fit
    d_eff <- ncol(se_sub) - ncol(design)
  }

  rowData(.data)[["dispersion_shrinked"]] <- as.numeric(disp)
  rowData(.data)[["dispersion_trended"]] <- as.numeric(trended)
  rowData(.data)[["dispersion_degrees_freedom"]] <-
    rep_len(as.numeric(d_eff), nrow(.data))
  .data <- attach_to_metadata(.data, fit, "estimateDisp")
  memorise_methods_used(.data, "edger")
}
