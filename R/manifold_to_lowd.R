# Author: Tyler J Burns
# Procedure: taking KNN of the original manifold, doing a low-d embedding,
#   and taking KNN of that.

#' @import Sconify
#' @importFrom magrittr "%>%"
#' @import Rtsne
#' @importFrom stats prcomp
#' @import tibble
#' @importFrom dplyr bind_cols
#' @importFrom stats predict
#' @import umap
NULL

#' @title Compare Neighborhoods
#' @description Calculates the sum of the union of each cell's KNN, one of
#' which comes from one method (eg original manifold) and the other of which
#' comes from another method (eg lowd embedding)
#' @param nn1 A matrix of cells by nearest neighborhood, with members being the
#' identify of a given nearest neighborhood. For example, an ID of 5 would
#' correspond to the fifth row of the original data matrix.
#' @param nn2 A matrix of cells by nearest neighborhood with the same format
#' as nn1, but the KNN were generated in a different way or on a different
#' piece or representation of the data.
#' @return A vector with each cell corresponding to its position in the
#' original data matrix. The values are the number of cells within the 2 KNN
#' that were in agreement, which in turn can be converted into a percent
#' fidelity
#' @export
CompareNeighborhoods <- function(nn1, nn2) {
    result <- sapply(1:nrow(nn1), function(i) {
        length(intersect(nn1[i,], nn2[i,]))
    })
    return(result)
}

#' @title Comparison Pipeline
#' @description Loop to generate KNN from original space, and lowd space (or
#' any other input), and determine how well the tSNE and PCA space approximate
#' the orignal manifold.
#' @param orig Tibble of cells by features
#' @param input_markers Vector of strings of names of surface markers of interest
#' @param lowd Tibble of dimr output
#' @param lowd_markers Vector of strings corresponding to the names of the lowd
#' embedding columns in your cells tibble.
#' @param k_titration Vector of values corresponding to the number of nearest
#' neighbors k to be tried. We recommend a range from very small (eg. 5) up to
#' half the dataset.
#' @return A tibble of cells by k values, where each value is a given cell's
#' local fidelity relative to the k value of interst.
#' @examples
#' library(dplyr)
#' k.titration <- c(10, 100)
#' tsne_names <- names(samusik_tsne)
#' ComparisonPipeline(samusik_cells, samusik_surface_markers, samusik_tsne, k.titration)
#' @export
ComparisonPipeline <- function(orig, lowd, input_markers, k_titration) {

    # Edge case for single k value
    if(length(k_titration) == 1) {
        nn_orig <- Fnn(cell.df = orig, input.markers = input_markers, k = k_titration)[[1]]
        nn_lowd <- Fnn(cell.df = lowd, input.markers = names(lowd), k = k_titration)[[1]]
        hl_compare <- CompareNeighborhoods(nn_orig, nn_lowd)/k_titration
        result <- as_tibble(hl_compare)
        names(result) <- k_titration
        return(result)
    }

    # Make the biggest value of K first
    k_titration <- sort(k_titration, decreasing = TRUE)

    # The KNN generation using the fnn command from Sconify
    nn_orig <- Fnn(cell.df = orig, input.markers = input_markers, k = k_titration[1])[[1]]
    nn_lowd <- Fnn(cell.df = lowd, input.markers = names(lowd), k = k_titration[1])[[1]]

    # Lowd compared to original manifold
    hl_compare <- CompareNeighborhoods(nn_orig, nn_lowd)/k_titration[1]
    hl_compare <- as_tibble(hl_compare)
    names(hl_compare) <- k_titration[1]
    # Return a tibble that has the percent average shared KNN for
    #   lowd versus original space

    # So there's no need to recompute KNN
    master_result <- lapply(k_titration[-1], function(i) {
        # A tracker
        message(i)
        nn_orig <- nn_orig[,1:i]
        nn_lowd <- nn_lowd[,1:i]
        hl_compare <- CompareNeighborhoods(nn_orig, nn_lowd)/i
        return(hl_compare)
    })

    # Turn this list into a tibble of cells by neighborhood size k
    names(master_result) <- k_titration[-1]
    master_result <- do.call(cbind, master_result) %>% as_tibble()

    # Create the final result and flip the columns back to original k titration
    master_result <- bind_cols(hl_compare, master_result) %>% .[,order(ncol(.):1)]
    return(master_result)
}

#' @title Run Principal Components Analysis (PCA)
#' @description Wrapper for 2-dimensional PCA
#' @param cells Tibble of cells by features
#' @param input Vector of markers to be considered as input for the PCA method
#' @return A tibble of cells by PC1 and PC2
#' @examples
#' RunPca(samusik_cells[1:1000,], samusik_surface_markers)
#' @export
RunPca <- function(cells, input) {
    pca <- prcomp(x = cells[,input])$x[,1:2] %>% as.tibble()
    return(pca)
}

#' @title Run t-SNE
#' @description Wrapper for 2-dimensional t-SNE
#' @param cells Tibble of cells by features
#' @param input Vector of markers to be considerd as input for the t-SNE method
#' @param perp The perplexity for t-SNE. Set to the CyTOF default of 30
#' @param to_pca Whether or not to PCA. Default = FALSE.
#' @return A tibble of cells by t-SNE1 and t-SNE2
#' @examples
#' RunTsne(samusik_cells[1:1000,], samusik_surface_markers)
#' @export
RunTsne <- function(cells, input = names(cells), perp = 30, to_pca = FALSE) {
    result <- Rtsne(X = cells[,input], perplexity = perp, verbose = TRUE, pca = to_pca, check_duplicates = FALSE)$Y %>% as.tibble()
    names(result) <- c("bh-SNE1", "bh-SNE2")
    return(result)
}

#' @title Run UMAP
#' @description Wrapper for 2-dimensional UMAP
#' @param cells Tibble of cells by features
#' @param input Vector of markers to be considered as input for UMAP method
#' @param nn Number of nearest neighbors to use in algorithm
#' @return A tibble of cells by UMAP1 and UMAP2
#' @examples RunUMAP(samusik_cells[1:1000,], samusik_surface_markers)
#' @export
RunUMAP <- function(cells, input = names(cells), nn = 15) {
    result <- umap::umap(cells[,input], n_neighbors = nn)$layout %>% as_tibble()
    names(result) <- c("umap1", "umap2")
    return(result)
}



