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
- **Optional arguments:** `distance` specifies the window size of QTLs to be tested. The default value is `distance=500000`, which includes SNPs within 500kb upstream and 500kb downstream to the CpG site. `b_vector` specifies the degree of shrinkage of each cell type. The default is `b_vector="0.2,0.2,0.2,0.2,0.2"` (the length must equal to the number of cell types).


If CTS data (in a small group of samples) is available, to incorporate the prior, please run the following 2 steps:

1. Create the prior data from CTS CpG-SNP associations, which have been obtained using the CTS data:
```
Rscript prepare_prior.r [file for priors] [file for cell type proportions] [path for outputs] \
Npair=1000000
```
- **file for priors:** The first five columns should be `probe`,`snp`,`chr`,`REF`,`ALT`. The remaining columns should be p-values and beta for each cell type with CTS data (not necessarily include all cell types in the file for cell type proportions). Column names are needed.
```
       probe       snp chr REF ALT      p_CD4T beta_CD4T      p_CD8T beta_CD8T
1 cg00045070 rs2479409   1   A   G 1.13095e-12   7.22779 2.95211e-14   7.37828
2 cg00345083 rs7517857   1   A   G 6.71566e-12   5.92522 3.45804e-10   5.24358
3 cg02890259  rs945417   1   G   C 3.13387e-13   7.97823 1.25687e-12   6.59110
4 cg02890259  rs945420   1   A   G 3.13387e-13   7.97823 1.25687e-12   6.59110
5 cg02890259  rs945421   1   C   T 3.13387e-13   7.97823 1.25687e-12   6.59110
6 cg02890259 rs6669935   1   G   C 1.05970e-09   7.48507 3.93825e-10   6.50018
...
```
- **file for cell type proportions:** This is the same input as that in `HBI.R`.
- **Optional arguments:** `Npair` specifies the total number of SNP-CpG pairs in the prior data, and the default value is the number of rows of the prior data. `p_adj_method` specifies the method used to adjust the prior data, and the default is `p_adj_method=bonferroni`. Other options can be `p_adj_method=fdr`,`p_adj_method=BY`.

2. Run the HBI algorithm with prior
```
Rscript HBI_cts_prior.R [file for phenotype/DNA methylation] [Phenotype/CpG name] [file for genotype] [file for cell type proportions] [file for covariates] [file for adjusted priors] [path for outputs] \
distance=500000
```
- **file for adjusted priors:** This is the output for the step 1 `prepare_prior.r`.
```
       probe       snp chr REF ALT      p_CD4T beta_CD4T      p_CD8T beta_CD8T
1 cg00045070 rs2479409   1   A   G 1.13095e-12   7.22779 2.95211e-14   7.37828
2 cg00345083 rs7517857   1   A   G 6.71566e-12   5.92522 3.45804e-10   5.24358
3 cg02890259  rs945417   1   G   C 3.13387e-13   7.97823 1.25687e-12   6.59110
4 cg02890259  rs945420   1   A   G 3.13387e-13   7.97823 1.25687e-12   6.59110
5 cg02890259  rs945421   1   C   T 3.13387e-13   7.97823 1.25687e-12   6.59110
6 cg02890259 rs6669935   1   G   C 1.05970e-09   7.48507 3.93825e-10   6.50018
   p_adj_CD8T beta_adj_CD8T  p_adj_CD4T beta_adj_CD4T
1 2.95211e-09      7.378280 1.13095e-07      7.227789
2 3.45804e-05      5.243399 6.71566e-07      5.925216
3 1.25687e-07      6.591099 3.13387e-08      7.978230
4 1.25687e-07      6.591099 3.13387e-08      7.978230
5 1.25687e-07      6.591099 3.13387e-08      7.978230
6 3.93825e-05      6.499924 1.05970e-04      7.484277

```
- Other inputs and optional arguments are the same as those in `HBI.R`.

## Credits
More information can be found in [our manuscript](https://www.biorxiv.org/content/10.1101/2024.02.01.578272v1).
