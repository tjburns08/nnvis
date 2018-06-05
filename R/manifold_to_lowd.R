# Author: Tyler J Burns
# Procedure: taking KNN of the original manifold, doing a low-d embedding,
#   and taking KNN of that.

#' @importFrom Sconify ProcessMultipleFiles Fnn
#' @importFrom magrittr "%>%"

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
#' @param cells Tibble of cells by features that contains also the lowd
#' embedding
#' @param input.markers The names of surface markers of interest
#' @param lowd.names String of vectors corresponding to the names of the lowd
#' embedding columns in your cells tibble.
#' @param k.titration Vector of values corresponding to the number of nearest
#' neighbors k to be tried. We recommend a range from very small (eg. 5) up to
#' half the dataset.
#' @return A tibble of cells by k values, where each value is a given cell's
#' local fidelity relative to the k value of interst.
#' @export
ComparisonPipeline <- function(cells, input.markers, lowd.names, k.titration) {
    master.result <- lapply(k.titration, function(i) {
        # A tracker
        message(i)

        # The KNN generation using the fnn command from Sconify
        nn.orig <- Fnn(cell.df = cells, input.markers = surface, k = i)[[1]]
        nn.lowd <- Fnn(cell.df = cells, input.markers = lowd.names, k = i)[[1]]

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

