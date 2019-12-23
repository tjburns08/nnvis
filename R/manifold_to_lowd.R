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
#' @param orig Tibble of cells by features that contains also the lowd
#' embedding
#' @param input.markers The names of surface markers of interest
#' @param lowd String of vectors corresponding to the names of the lowd
#' embedding columns in your cells tibble.
#' @param k.titration Vector of values corresponding to the number of nearest
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
ComparisonPipeline <- function(orig, input.markers, lowd, k.titration) {
    master.result <- lapply(k.titration, function(i) {
        # A tracker
        message(i)

        # The KNN generation using the fnn command from Sconify
        nn.orig <- Fnn(cell.df = orig, input.markers = input.markers, k = i)[[1]]
        nn.lowd <- Fnn(cell.df = lowd, input.markers = names(lowd), k = i)[[1]]

        # Lowd compared to original manifold
        hl.compare <- CompareNeighborhoods(nn.orig, nn.lowd)/i

        # Return a tibble that has the percent average shared KNN for
        #   lowd versus original space
        return(hl.compare)
    })

    # Turn this list into a tibble of cells by neighborhood size k
    names(master.result) <- k.titration
    master.result <- do.call(cbind, master.result) %>% as.tibble()
    return(master.result)
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
    result <- Rtsne(X = cells[,input], perplexity = perp, verbose = TRUE, pca = to_pca)$Y %>% as.tibble()
    names(result) <- c("bh-SNE1", "bh-SNE2")
    return(result)
}

#' @title Run UMAP
#' @description Wrapper for 2-dimensional UMAP
#' @param cells Tibble of cells by features
#' @param input Vector of markers to be considered as input for UMAP method
#' @return A tibble of cells by UMAP1 and UMAP2
#' @examples RunUMAP(samusik_cells[1:1000,], samusik_surface_markers)
#' @export
RunUMAP <- function(cells, input = names(cells)) {
    result <- umap(cells[,input])$layout
    names(result) <- c("umap1", "umap2")
    return(result)
}



