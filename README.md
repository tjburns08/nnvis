# nnvis
Make KNN-based identity comparisons between different manifolds (eg.  original space vs t-SNE space).

## Description
This R package allows one to ask and answer the question "how well is this t-SNE map (or other low-dimensional embedding) 
representing my data." This is done by simply taking the K-nearest neighbors of the original manifold, and the k-nearest 
neighbors of the lower dimensional embedding, and comparing them for each cell. These values can be averaged across a wide 
range of values k to assess global similarity, or visualized directly on a low-dimensional embedding to assess local similarity. 

## How to use
Get your data in a data frame or tibble format of cells x features. Run your low-dimensional embedding of interest, name
the parameters, and merge them with your original data. Then run the ComparisonPipeline function. This will output a tibble
of cells by features in the same order as the original data, with each element being the KNN similarity between original 
space and embedded space for that particular cell, as a percentage. Now you can color a low-dimensional embedding by its 
own efficacy!
