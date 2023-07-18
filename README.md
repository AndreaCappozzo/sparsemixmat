
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sparsemixmat

<!-- badges: start -->
<!-- badges: end -->

The goal of sparsemixmat is to provide a framework for matrix-variate
data sparse model-based clustering. The approach assumes that the matrix
mixture parameters are sparse and have different degree of sparsity
across clusters, allowing to induce parsimony in a flexible manner.
Estimation of the model relies on the maximization of a penalized
likelihood, with specifically tailored group and graphical lasso
penalties.

This repository is associated with the paper Cappozzo, Casa, Fop (2023+)
Sparse model-based clustering of three-way data via lasso-type
penalties. FIXME ADD ARXIV LINK

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("AndreaCappozzo/sparsemixmat")
```
