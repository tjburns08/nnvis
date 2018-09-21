#' Publically available dataset from Nikolay Samusik of the Nolan Lab,
#' consisting of mouse bone marrow. Original data is 86,864 cells, but here
#' it is subsampled down to 10,000. The data were asinh transformed. The
#' Flow Repository ID FR-FCM-ZZPH.
#' @format a tibble of 10000 cells by 54 features
"samusik_cells"

#' Selected surface markers from the dqvis_cells dataset to be used for dim
#' reduction analysis
#' @format a vector of 38 features
"samusik_surface_markers"

#' The output of t-SNE results from the dqvis_cells dataset with the
#' dqvis_surface_markers used as input.
#' @format a tibble of 10000 cells by 2 features
"samusik_tsne"
