#' Plot reduced dimensions
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description plot_reduced_dimension() takes a `SummarizedExperiment` that has
#' been processed with `reduce_dimensions()` and returns a ggplot of the selected
#' reduced dimensions, optionally coloured by a sample covariate. For PCA, axis
#' labels include the percentage of variance explained.
#'
#' @importFrom rlang enquo sym
#' @importFrom dplyr case_when
#' @importFrom stringr str_to_lower str_detect
#' @importFrom S4Vectors metadata
#' @importFrom ggplot2 ggplot aes geom_point labs
#'
#' @name plot_reduced_dimension
#'
#' @param .data A `SummarizedExperiment` object
#' @param .color A symbol. Column in `colData` used for colour coding
#' @param method A character string. Dimensionality reduction method to plot
#'   (e.g. `"PCA"`, `"MDS"`, `"tSNE"`, `"UMAP"`)
#' @param dims Integer vector of length 2. Which dimensions to plot (default `1:2`)
#'
#' @return A ggplot object
#'
#' @examples
#' ## Load airway dataset for examples
#' data('airway', package = 'airway')
#'   # Ensure a 'condition' column exists for examples expecting it
#'
#'     SummarizedExperiment::colData(airway)$condition <- SummarizedExperiment::colData(airway)$dex
#'
#'
#' counts.PCA =
#'  airway |>
#'  identify_abundant() |>
#'  reduce_dimensions(assay = "counts", method="PCA", .dims = 3)
#'
#'  plot_reduced_dimension(counts.PCA, .color = condition, method = "PCA", dims = 1:2)
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' @docType methods
#' @rdname plot_reduced_dimension-methods
#' @export
#'
setGeneric("plot_reduced_dimension", function(.data,
                                              .color,
                                              method = "PCA",
                                              dims = 1:2)
  standardGeneric("plot_reduced_dimension"))

.plot_reduced_dimension_se <- function(.data,
                                       .color,
                                       method = "PCA",
                                       dims = 1:2) {
  .color <- enquo(.color)

  if (length(dims) != 2)
    stop("tidybulk says: `dims` must be an integer vector of length 2.", call. = FALSE)

  dim_type <- case_when(
    str_to_lower(method) == str_to_lower("PCA") ~ "pca",
    str_to_lower(method) == str_to_lower("MDS") ~ "mds",
    str_to_lower(method) == str_to_lower("tSNE") ~ "tsne",
    str_to_lower(method) == str_to_lower("UMAP") ~ "umap",
    TRUE ~ NA_character_
  )

  if (is.na(dim_type))
    stop(
      "tidybulk says: method must be one of \"PCA\", \"MDS\", \"tSNE\", or \"UMAP\".",
      call. = FALSE
    )

  coord_naming <- case_when(
    dim_type == "pca" ~ "PC",
    dim_type == "mds" ~ "Dim",
    dim_type == "tsne" ~ "tSNE",
    dim_type == "umap" ~ "UMAP"
  )

  df <- .data |> pivot_sample()

  coord_names <- names(df)[str_detect(names(df), paste0("^", coord_naming, "\\d+"))]
  if (length(coord_names) < max(dims))
    stop(
      "tidybulk says: could not find enough \"",
      coord_naming,
      "*\" columns for the requested dims. Did you run reduce_dimensions() with method=\"",
      method,
      "\"?",
      call. = FALSE
    )
  coord_names <- coord_names[dims]

  p <- df |>
    ggplot(aes(!!sym(coord_names[1]), !!sym(coord_names[2]), color = !!.color)) +
    geom_point()

  if (dim_type == "pca") {
    pca_res <- metadata(.data)$tidybulk$PCA
    vars <- pca_res[[1]]^2 / sum(pca_res[[1]]^2)

    p <- p +
      labs(
        x = paste0(coord_names[1], " (", round(vars[dims[1]] * 100), "% variance explained)"),
        y = paste0(coord_names[2], " (", round(vars[dims[2]] * 100), "% variance explained)")
      )
  }

  p
}

#' plot_reduced_dimension
#'
#' @docType methods
#' @rdname plot_reduced_dimension-methods
#'
#' @return A ggplot object
#'
setMethod("plot_reduced_dimension",
          "SummarizedExperiment",
          .plot_reduced_dimension_se)

#' plot_reduced_dimension
#'
#' @docType methods
#' @rdname plot_reduced_dimension-methods
#'
#' @return A ggplot object
#'
setMethod("plot_reduced_dimension",
          "RangedSummarizedExperiment",
          .plot_reduced_dimension_se)
