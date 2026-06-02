#' Plot reduced dimensions
#'
#' @param SE A 'SummarizedExperiment'
#' @param VAR A Symbol. What variable in colData to be used for color coding
#' @param method A character string. IndicatingMethod of dimensionality reduction to be ploted 
#' @param dims Vector of dimensions to be plotted
#'
#' @returns A ggplot object
#' @export
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
#'  plot_reduced_dim(counts.PCA, condition, "PCA", 1:2)
plot_reduced_dim <- function(SE, VAR, method = "PCA", dims = 1:2){
  varname <- rlang::ensym(VAR)
  
  dim_type <- dplyr::case_when(
    stringr::str_to_lower(method) == stringr::str_to_lower("PCA") ~ "pca",
    stringr::str_to_lower(method) == stringr::str_to_lower("MDS") ~ "mds",
    stringr::str_to_lower(method) == stringr::str_to_lower("tSNE") ~ "tsne",
    stringr::str_to_lower(method) == stringr::str_to_lower("UMAP") ~ "umap")
  
  coord_naming <- dplyr::case_when(
    dim_type == "pca" ~ "PC",
    dim_type == "mds" ~ "Dim",
    dim_type == "tsne" ~ "tSNE",
    dim_type == "umap" ~ "UMAP")
  
  df <- SE %>% 
    pivot_sample()
  
  coord_names <- names(df)[stringr::str_detect(names(df), paste0("^", coord_naming, "\\d+"))]
  coord_names <- coord_names[dims]
  
  
  p <- SE |>
    pivot_sample() |>
    ggplot2::ggplot(ggplot2::aes(.data[[coord_names[1]]], .data[[coord_names[2]]], color = !!varname))+
    ggplot2::geom_point()
  
  if(dim_type == "pca"){
    pca_res <- S4Vectors::metadata(SE)$tidybulk$PCA
    vars <- pca_res[[1]]^2 / sum(pca_res[[1]]^2)
    
    p <- p +
      ggplot2::labs(x = paste0("PC1 (", round(vars[1] * 100), "% variance explained)"),
                    y = paste0("PC2 (", round(vars[2] * 100), "% variance explained)"))
  }
  return(p)
}