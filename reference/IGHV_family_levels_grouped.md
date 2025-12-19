# Human IGHV families in ascending numeric order with IGHV6 and 7 grouped

A vector that can be used to select and arrange IGHV families while
IGHV6 and 7 are grouped as they only contain one gene each. This
requires that IGHV_family is set to IGHV6/7 for those families, which is
done automatically by TS_plot_IG_usage(compare_what =
"IGHV_family_grouped") Data obtained from IMGT

## Usage

``` r
as.factor(string, levels = TS::IGHV_family_levels_grouped)
```

## Format

A character vector of length 6
