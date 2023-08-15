# DepMap
RNA Expression and Gene Dependency Analysis to Select Cell Lines for CRISPR knock out
Code is written to determine RNA expression of genes and Gene dependency values from DepMap (Nusinow et al 2020)
Multiple data transformation steps for data downloaded from depmap were needed to run the code
Most recent depmap data were downloaded from https://depmap.org/portal/download/all/?releasename=DepMap+Public+23Q2&filename=Model.csv
Cut off for RNA expression is Log2TMP+1>3 but can be chagned in the For loop
Reference for dataset used: Nusinow, D. P., J. Szpyt, M. Ghandi, C. M. Rose, E. R. McDonald, 3rd, M. Kalocsay, J. Jane-Valbuena, E. Gelfand, D. K. Schweppe, M. Jedrychowski, J. Golji, D. A. Porter, T. Rejtar, Y. K. Wang, G. V. Kryukov, F. Stegmeier, B. K. Erickson, L. A. Garraway, W. R. Sellers and S. P. Gygi (2020). "Quantitative Proteomics of the Cancer Cell Line Encyclopedia." Cell 180(2): 387-402 e316.
