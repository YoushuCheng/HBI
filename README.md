# HBI
A hierarchical Bayesian interaction model to estimate cell-type-specific methylation quantitative trait loci incorporating priors from cell-sorted bisulfite sequencing data

## Tutorial
Please run Rscript HBI.R [file for DNA methylation] [file for genotype] [file for cell type proportions] [file for covariates] [file for priors] [path for outputs]

- Expression file should not contain col.names and row.names

- Batch index file should not contain col.names and row.names

- Biological group index file should not contain col.names and row.names

The output imputed and de-noised expression matrix is in the same size of the input expression matrix.

## Credits
More information can be found in [our manuscript](https://www.biorxiv.org/content/10.1101/2024.02.01.578272v1).
