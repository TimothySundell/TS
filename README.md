
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TS - A collection of QoL functions <a href="https://github.com/TimothySundell/TS"><img src="man/Figures/TS_logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![](https://img.shields.io/github/languages/code-size/TimothySundell/TS.svg)](https://github.com/TimothySundell/TS)
[![](https://img.shields.io/badge/doi-10.1093/bfgp/elac044-green.svg)](https://doi.org/10.1093/bfgp/elac044)
<!-- badges: end -->

TS was written to improve the quality of life for its creator, maybe it
can help you as well?

It contains a collection of functions that may be relevant people
analysing single-cell RNA-sequencing data, and analysing
Ig(BCR)/TCR-sequences.

**Please** give our method article a look at the **DOI** button above.
It shows the importance of removing BCR/TCR genes prior to performing
unsupervised clustering and downstream analyses of scRNA-seq data
containing lymphocytes.

You can find a reference webpage
[*here*](https://timothysundell.github.io/TS/), or under the *About*
section in the top right corner of this webpage.

------------------------------------------------------------------------

## Installation

To install the latest version of the TS package:

``` r
devtools::install_github("TimothySundell/TS")
```

------------------------------------------------------------------------

## Comments from the author

My plan is to keep this repository updated as much as possible. If
something is not working, please post in the *Issues* section.

As of now, the package does not pass *R CMD CHECK* as the tidyverse way
of coding leaves me a lot of unlinked global variables to reference.
Will maybe be fixed in the future.

------------------------------------------------------------------------

## The TS collection of functions can be divided into 4 main categories:

- Functions for plotting figures
- Functions for analysing sc-V(D)J-sequencing data
- Functions for analysing scRNA-seq data
- General data-wrangling

## Repo activity

<figure>
<img
src="https://repobeats.axiom.co/api/embed/6ceecc56e6072d87a2a5cba22d6ddb12968a964d.svg"
title="Repobeats analytics image" alt="Repo activity" />
<figcaption aria-hidden="true">Repo activity</figcaption>
</figure>
