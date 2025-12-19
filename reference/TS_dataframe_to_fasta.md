# Converts a dataframe and saves it in Fasta format

Requires a dataframe with two columns:

## Usage

``` r
TS_dataframe_to_fasta(
  dataframe,
  filename = gsub(pattern = " ", replacement = "_", x = gsub(pattern = ":", replacement =
    "-", x = paste0(Sys.time(), ".fasta")))
)
```

## Arguments

- dataframe:

  A 2-column dataframe where:

  - Column 1 contains cell IDs

  - Column 2 contains the full nucleotide sequence of the contig
    (concatenated fwr and cdr sequences in nt-form)

- filename:

  Filename for the output FILE.fasta. Defaults to 'system date and
  time.fasta'
