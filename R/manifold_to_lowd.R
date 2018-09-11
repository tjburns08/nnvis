# Author: Tyler J Burns
# Procedure: taking KNN of the original manifold, doing a low-d embedding,
#   and taking KNN of that.

#' @import Sconify
#' @importFrom magrittr "%>%"
#' @import Rtsne
#' @import umap
#' @import keras
#' @importFrom stats prcomp
#' @import tibble
#' @importFrom dplyr bind_cols

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
#' @examples
#' library(dplyr)
#' k.titration <- c(10, 100)
#' tsne_names <- names(dqvis_tsne)
#' ComparisonPipeline(dqvis_cells, dqvis_tsne, tsne_names, k.titration)
#' @export
ComparisonPipeline <- function(orig, lowd, k.titration) {
    master.result <- lapply(k.titration, function(i) {
        # A tracker
        message(i)

        # The KNN generation using the fnn command from Sconify
        nn.orig <- Fnn(cell.df = orig, input.markers = names(orig), k = i)[[1]]
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
#' RunPca(dqvis_cells[1:1000,], dqvis_surface_markers)
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
#' @return A tibble of cells by t-SNE1 and t-SNE2
#' @examples
#' RunTsne(dqvis_cells[1:1000,], dqvis_surface_markers)
#' @export
RunTsne <- function(cells, input = names(cells), perp = 30) {
    result <- Rtsne(X = cells[,input], perplexity = perp, verbose = TRUE)$Y %>% as.tibble()
    names(result) <- c("bh-SNE1", "bh-SNE2")
    return(result)
}

#' @title Runs UMAP on one's cells
#' @description Wrapper for 2-dimensional UMAP
#' @param cells TIbble of cells by features
#' @param input Vector of markers to be considered as input for the UMAP method
#' @return A tibble of cells by UMAP1 and UMAP2
#' @examples
#' RunUmap(dqvis_cells[1:1000,])
#' @export
RunUmap <- function(cells, input = names(cells)) {
    result <- umap(as.matrix(cells[,input]))$layout
    return(result)
}

#' @title Run variational autoencoder
#' @description Wrapper for simple variational autoencoder through keras
#' @param cells: Tibble of cells by features
#' @param input Vector of markers to be considered as input for the VAE
#' @return A tibble of cells by enc1 and enc2
#' @examples
#' RunVae(cells[1:1000,])
#' @export
RunVae <- function(cells, input = names(cells), epochs = 50L) {
    # Code below re-purposed from:
    # https://keras.rstudio.com/articles/examples/variational_autoencoder.html
    # Parameters --------------------------------------------------------------
    K <- keras::backend()
    batch_size <- 100L
    original_dim <- length(input) %>% as.integer()
    latent_dim <- 2L
    intermediate_dim <- 8L
    epochs <- epochs %>% as.integer()
    epsilon_std <- 1.0

    # Model definition -------------------------------------------------------

    x <- layer_input(shape = c(original_dim))
    h <- layer_dense(x, intermediate_dim, activation = "relu")
    z_mean <- layer_dense(h, latent_dim)
    z_log_var <- layer_dense(h, latent_dim)

    sampling <- function(arg){
        z_mean <- arg[, 1:(latent_dim)]
        z_log_var <- arg[, (latent_dim + 1):(2 * latent_dim)]

        epsilon <- k_random_normal(
            shape = c(k_shape(z_mean)[[1]]),
            mean=0.,
            stddev=epsilon_std
        )

        z_mean + k_exp(z_log_var/2)*epsilon
    }

    # note that "output_shape" isn't necessary with the TensorFlow backend
    z <- layer_concatenate(list(z_mean, z_log_var)) %>%
        layer_lambda(sampling)

    # we instantiate these layers separately so as to reuse them later
    decoder_h <- layer_dense(units = intermediate_dim, activation = "relu")
    decoder_mean <- layer_dense(units = original_dim, activation = "sigmoid")
    h_decoded <- decoder_h(z)
    x_decoded_mean <- decoder_mean(h_decoded)

    # end-to-end autoencoder
    vae <- keras_model(x, x_decoded_mean)

    # encoder, from inputs to latent space
    encoder <- keras_model(x, z_mean)

    # generator, from latent space to reconstructed inputs
    decoder_input <- layer_input(shape = latent_dim)
    h_decoded_2 <- decoder_h(decoder_input)
    x_decoded_mean_2 <- decoder_mean(h_decoded_2)
    generator <- keras_model(decoder_input, x_decoded_mean_2)


    vae_loss <- function(x, x_decoded_mean){
        xent_loss <- (original_dim/1.0)*loss_binary_crossentropy(x, x_decoded_mean)
        kl_loss <- -0.5*k_mean(1 + z_log_var - k_square(z_mean) - k_exp(z_log_var), axis = -1L)
        xent_loss + kl_loss
    }

    vae %>% compile(optimizer = "rmsprop", loss = vae_loss)

    # Defining the cells
    cells <- cells[,input]

    # MinMax normalization
    cells <- apply(cells, 2, function(i) {
        normalized = (x-min(i))/(max(i)-min(i))
    })

    # Training data and test data will be the same right here
    x_train <- cells
    x_test <- cells

    # Initialize model
    #InitializeVae()
    # From above, start with fresh model
    vae <- keras_model(x, x_decoded_mean)
    vae %>% compile(optimizer = "rmsprop", loss = vae_loss)

    # Autoencoder model
    vae %>% fit(
        x_train, x_train,
        shuffle = TRUE,
        epochs = epochs,
        batch_size = batch_size,
        validation_data = list(x_test, x_test)
    )

    # Generate latent dimensions for your data
    x_test_encoded <- predict(encoder, x_test, batch_size = batch_size)
    colnames(x_test_encoded) <- c(paste("enc1", i, sep = "_"), paste("enc2", i, sep = "_"))
    return(list(result = x_test_encoded, encoder = encoder, batch_size = batch_size))
}

