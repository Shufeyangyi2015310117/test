# SpatialAnno
SpatialAnno: Probabilistic cell/domain-type assignment of spatial transcriptomics data with SpatialAnno

SpatialAnno is a package for annotation on spatial transcriptomics datasets developed by Jin Liu's lab. It has the capability to effectively leverage a large number of non-marker genes as well as “qualitative” information about marker genes without using a reference dataset. Uniquely, SpatialAnno estimates low-dimensional embeddings for a large number of non-marker genes via a factor model while promoting spatial smoothness among neighboring spots via a Potts model

# Installation

To install the packages "SpatialAnno", firstly, install the 'devtools' package. Besides, "SpatialAnno" depends on the 'Rcpp' and 'RcppArmadillo' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

`install.packages("devtools")`

`library(devtools)`

`install_github("Shufeyangyi2015310117/SpatialAnno")`

It requires a few minutes on a "normal" desktop computer.  


# Demonstration

For an example of typical SpatialAnno usage, please see our [Package vignette](https://shufeyangyi2015310117.github.io/SpatialAnno/index.html) for a demonstration and overview of the functions included in SpatialAnno.

# Usage
For usage examples and guided walkthroughs, check the `vignettes` directory of the repo. 

* [Mouse olfactory bulb data analysis](https://shufeyangyi2015310117.github.io/SpatialAnno/articles/MOB.html)
* [Toy examples for Annotation](https://shufeyangyi2015310117.github.io/SpatialAnno/articles/SpatialAnno.html)
* [human dorsolateral prefrontal cortex data analysis](https://shufeyangyi2015310117.github.io/SpatialAnno/articles/brain.html)

# Analysis code

The analysis code of SpatialAnno are accessable on [code website](https://github.com/Shufeyangyi2015310117/SpatialAnno_Analysis)

We leave emply the folders containing inputs and intermediate outputs due to the size limit. The whole files can be downloaded from [zenodo](https://doi.org/10.5281/zenodo.7413083)
