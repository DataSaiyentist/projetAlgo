# projetAlgo - Projet Algorithmique (M2 Data Science 2022/2023)

Le projet est d’implémenter FOCuS en R/C++ et de comparer ses performances (algorithmiques et statistiques) à celle de l’approche séquentielle de Page.
</br>
Ces méthodes permettent de détecter des ruptures dans une séquence (en particulier des séries temporelles).

## Installation du package

Installez le package sur R de la manière suivante :

```{r gitInstall, eval = FALSE}
devtools::install_github("DataSaiyentist/projetAlgo")
library(projetAlgo)
```

## References

- Romano, Gaetano and Eckley, Idris and Fearnhead, Paul and Rigaill, Guillem (2021) Fast Online Changepoint Detection via Functional Pruning CUSUM statistics, arXiv, [10.48550/ARXIV.2110.08205](https://arxiv.org/abs/2110.08205)
- [https://github.com/gtromano/FOCuS](https://github.com/gtromano/FOCuS)

## License

Copyright © 2021 [Data Saiyentist](https://github.com/DataSaiyentist). <br />
This project is [GNU General Public License v3.0]() licensed.