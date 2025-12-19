# Compute pairwise Fisher tests for gene usage

Given a data.frame of gene counts per group, performs Fisher's Exact
Test for all pairwise group combinations for each gene.

## Usage

``` r
TS_compute_pairwise_fisher(df, gene_var, group_var, p_adj_method = "hochberg")
```

## Arguments

- df:

  A data.frame or tibble containing at least two columns: one for genes
  and one for groups.

- gene_var:

  String; name of the column in `df` that identifies the gene (e.g.
  "gene_list").

- group_var:

  String; name of the column in `df` that identifies the group (e.g.
  "named_clusters").

- p_adj_method:

  String; method for multiple-testing correction passed to `p.adjust`
  (default: "hochberg").

## Value

A tibble with columns: gene, group1, group2, p_value, p_adj, p_sig
