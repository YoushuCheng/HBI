# HBI
A hierarchical Bayesian interaction model to estimate cell-type-specific methylation quantitative trait loci incorporating priors from cell-sorted bisulfite sequencing data

## Tutorial
If only bulk data is available, please run 
```
Rscript HBI.R [file for phenotype/DNA methylation] [Phenotype/CpG name] [file for genotype] [file for cell type proportions] [file for covariates] [path for outputs] \
distance=500000
```
- **File for phenotype/DNA methylation:** Each row represent a phenotype (CpG site in DNA methylations), and the first three columns are `probe`,`BP`,`CHR`, the remaining columns represent individuals. Column names are needed.
```
        probe       BP CHR        id_1       id_2        id_3       id_4
1  cg08730728 37252593  22 -5.48097920 -5.7615145 -5.84497641 -5.7453667
2  cg25294651 42919443  22  2.92992156  2.7553017  2.30812230  2.3923174
3  cg15283028 41601160  22 -6.07666920 -5.9796267 -6.37562358 -6.2921390
4  cg21855135 42353933  22 -2.55987152 -2.1778738 -2.75530165 -2.6189098
5  cg04371950 45602057  22 -0.07503705  0.1965098 -0.09236402 -0.5789538
6  cg03070533 35653394  22 -6.30378075 -6.4641081 -6.42553520 -6.3879443
...
```
- **Phenotype/CpG name:** The CpG to be tested should be specified, for example, `cg08730728`.

If CTS data (in a small group of samples) is available, to incorporate the prior, please run the following 2 steps:

1. Create the prior data from CTS CpG-SNP associations, which have been obtained using the CTS data:
```
Rscript prepare_prior.r [file for priors] [file for cell type proportions] [path for outputs] \
Npair=1000000
```

2. Run the HBI algorithm with prior
```
Rscript HBI_cts_prior.R [file for phenotype/DNA methylation] [Phenotype/CpG name] [file for genotype] [file for cell type proportions] [file for covariates] [file for priors] [path for outputs] \
distance=500000
```

- Expression file should not contain col.names and row.names

- Batch index file should not contain col.names and row.names

- Biological group index file should not contain col.names and row.names

The output imputed and de-noised expression matrix is in the same size of the input expression matrix.

## Credits
More information can be found in [our manuscript](https://www.biorxiv.org/content/10.1101/2024.02.01.578272v1).
