# SpatialAnno
SpatialAnno: Probabilistic cell/domain-type assignment of spatial transcriptomics data with SpatialAnno

SpatialAnno is a package for annotation on spatial transcriptomics datasets developed by Jin Liu's lab. It has the capability to effectively leverage a large number of non-marker genes as well as “qualitative” information about marker genes without using a reference dataset. Uniquely, SpatialAnno estimates low-dimensional embeddings for a large number of non-marker genes via a factor model while promoting spatial smoothness among neighboring spots via a Potts model

# Installation

To install the packages "SpatialAnno", firstly, install the 'devtools' package. Besides, "SpatialAnno" depends on the 'Rcpp' and 'RcppArmadillo' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

`install.packages("devtools")`

`library(devtools)`

`install_github("Shufeyangyi2015310117/SpatialAnno")`

It requires a few minutes on a "normal" desktop computer.  


# Issues

For the issues met in the usage of SpatialAnno, please see our [Issues website](https://github.com/Shufeyangyi2015310117/SpatialAnno/issues) for help. You can also report a bug in that section. we will solve it as soon as possible.
