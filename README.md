# nnvis
Make KNN-based identity comparisons between different manifolds (eg.  original space vs t-SNE space).

## Description
This R package allows one to ask and answer the question "how well is this t-SNE map (or other low-dimensional embedding) 
representing my data." This is done by simply taking the K-nearest neighbors of the original manifold, and the k-nearest 
neighbors of the lower dimensional embedding, and comparing them for each cell. These values can be averaged across a wide 
range of values k to assess global similarity, or visualized directly on a low-dimensional embedding to assess local similarity. 

## Installation
```
# From GitHub
library(devtools)
devtools::install_github("tjburns08/nnvis")
```

## How to use
Start with the example data provided (Samusik_01.fcs). See the vignette for
code accordingly. 


