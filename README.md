# HBI
A hierarchical Bayesian interaction model to estimate cell-type-specific methylation quantitative trait loci incorporating priors from cell-sorted bisulfite sequencing data

## Tutorial
If only bulk data is available, please run 
```
Rscript HBI.R [file for phenotype/DNA methylation] [Phenotype/CpG name] [file for genotype] [file for cell type proportions] [file for covariates] [path for outputs] \
distance=500000
```
- **file for phenotype/DNA methylation:** Each row represent a phenotype (CpG site in DNA methylations), and the first three columns should be `probe`,`BP`,`CHR`, the remaining columns represent individuals. Column names are needed.
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
- **file for genotype:** Genotype file in VCF format.
- **file for cell type proportions:** The first column is `IID`, each of the remaining columns represents one cell type. Column names are needed.
```
   IID      CD8T       CD4T       Mono        NK      Bcell
1 id_1 0.2549020 0.24509804 0.08823529 0.2352941 0.17647059
2 id_2 0.2626263 0.45454545 0.02020202 0.1010101 0.16161616
3 id_3 0.4019608 0.05882353 0.06862745 0.2843137 0.18627451
4 id_4 0.3861386 0.26732673 0.04950495 0.1386139 0.15841584
5 id_5 0.2765957 0.25531915 0.15957447 0.2127660 0.09574468
6 id_6 0.3700000 0.27000000 0.09000000 0.1100000 0.16000000
...
```
- **file for covariates:** The first column is `IID`, each of the remaining columns represents one covariate (dummy variables should have be created for categorical variables). Column names are needed.
```
   IID         PC1         PC2         PC3         PC4          PC5 AGEATVIS
1 id_1 -0.07332110 -0.00368784 -0.01427090 -0.01352250 -0.015924600    63.36
2 id_2 -0.00802792  0.00339643  0.01105780  0.02708000 -0.006274620    24.57
3 id_3 -0.06502260  0.00238952 -0.00620466 -0.00135235 -0.002030430    48.07
4 id_4 -0.07117030  0.00157266 -0.00716070 -0.01088720  0.011414200    41.76
5 id_5  0.04568470  0.06306540 -0.01716520 -0.01318330 -0.000486778    50.22
6 id_6 -0.00226559  0.00160982  0.03022440  0.03769350 -0.003312670    47.74
...
```
- **path for outputs:** For example, `/Mypath/result.txt`.

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
