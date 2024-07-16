# HBI
A hierarchical Bayesian interaction model to estimate cell-type-specific methylation quantitative trait loci incorporating priors from cell-sorted bisulfite sequencing data

Please run Rscript HBI.R [expression file] [batch index file] [biological group index file] [output file]

- Expression file should not contain col.names and row.names

- Batch index file should not contain col.names and row.names

- Biological group index file should not contain col.names and row.names

The output imputed and de-noised expression matrix is in the same size of the input expression matrix.

More information can be found in [our manuscript](https://www.biorxiv.org/content/10.1101/2024.02.01.578272v1).
