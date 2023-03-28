# Projet d'Algorithmique (M2 Data Science 2022/2023)

L'objectif est d’implémenter des méthodes permettant de détecter des ruptures dans une séquence (en particulier des séries temporelles) comme CUSUM ou Page-CUSUM.

## Exigences pour l'installation du package

L'installation nécessite les packages suivants :

- `devtools`
- `Rcpp (>= 1.0.5)`

## Installation du package

Installez le package sur R de la manière suivante :

```r
devtools::install_github("DataSaiyentist/projetAlgo")
library(projetAlgo)
```

## Références

- Romano, Gaetano and Eckley, Idris and Fearnhead, Paul and Rigaill, Guillem (2021) Fast Online Changepoint Detection via Functional Pruning CUSUM statistics, arXiv, [10.48550/ARXIV.2110.08205](https://arxiv.org/abs/2110.08205)
- [https://github.com/gtromano/FOCuS](https://github.com/gtromano/FOCuS)

## License

Copyright © 2023 [Data Saiyentist](https://github.com/DataSaiyentist). <br />
Ce projet est sous license [GNU General Public License v3.0](https://github.com/DataSaiyentist/projetAlgo/blob/main/LICENSE).
