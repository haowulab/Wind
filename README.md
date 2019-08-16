# Wind: weighted indexes for clustering evaluation

`Wind` is an R package designed to compute weighted normalized mutual information (wNMI) and weighted Rand index (wRI) to evaluate the clustering results by comparing a clustering output with a reference which has a hierarchical structure.


## Installation
---------------

Run following commands in R :

```
library(devtools)
install_github("haowulab/Wind", build_opts = c("--no-resave-data"))
```

The second line might take a little time and install some extra packages for building the package vignette.

For detailed information about the package, please refer to the package vignette: 

```
library(Wind)
vignette("Wind")
```

