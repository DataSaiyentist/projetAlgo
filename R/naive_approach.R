######################################################################
###                                                                ###
###  GPL-3 License                                                 ###
###  Copyright (c) 2023 DataSaiyentist, Zaari1998, Ibtissem Rebai  ###
###                                                                ###
######################################################################

#' Test de Wilcoxon
#'
#' @description Cette fonction retourne la version corrigée de la p-value du test de Wilcoxon (ou somme des rangs) pour deux séquences non appariées.
#'
#' @details Sous l'hypothèse nulle, on suppose que les deux séquences proviennent de la même distribution et que \code{Mediane(Y-X) = 0}.
#'
#' @param x vecteur (série temporelle) à valeurs réelles
#' @param y vecteur (série temporelle) à valeurs réelles que l'on souhaite comparer à x
#' @return \code{p.value} la p-value du test de Wilcoxon
#' @seealso https://github.com/obouaziz/robusTest (Inspiration pour l'implémentation)
wilcoxtest  <- function(x, y) {

  n <- length(x)
  m <- length(y)
  ties <- x%in%y

  # Au cas où on aurait des valeurs égales entre x et y
  if (sum(ties)!=0) {
    # On ajoute une petite valeure aléatoire entre -1e-5 et 1e-5
    # pour permettre de mieux distinguer les deux séquences
    x[ties] <- x[ties] + runif(sum(ties), -0.00001, 0.00001)
  }

  # Calcul de la statistique de test
  R <- array(0, dim=c(n,m))
  for (i in 1:n) {
    for (j in 1:m) {
      R[i,j] <- (y[j]>x[i])
    }
  }

  H <- apply(R, 1, mean)
  G <- apply(R, 2, mean)
  V <- var(H)/n + var(G)/m
  Tn <- (mean(R) - 0.5)/sqrt(V)

  # Calcul de la p.value (On rejette Ho, si |Tn| > seuil)
  Pval <- 2*(1 - pnorm(abs(Tn)))

  return(list(p.value = Pval))
}

######################################################################
######################################################################

#' Première approche naïve
#'
#' @description Cette fonction teste sur chaque point (votre choix) le test de Wilcoxon ou de Kologorov-Smirnov.
#'
#' @details Elle retourne les indices des points pour lesquelles la p-value est plus petite q'un certain seuil.
#' Attention, le test de Kolmogorov-Smirnov n'a pas été implémenté à la main (\code{ks.test} de la librairie \code{"stats"}).
#'
#' @param x vecteur (série temporelle) à valeurs réelles
#' @param threshold seuil pour décider d'un changement de tendance
#' @param method statistique à utiliser parmi "wilcox" et "KS"
#' @return \code{index_break} les indices associées aux ruptures
naive_test <- function(x, threshold = 0.05, method = c("wilcox", "KS")) {

  n <- length(x)

  if (method=="wilcox") { test_func <- wilcoxtest }
  if (method=="KS") { test_func <- ks.test }

  pvalues <- NULL
  # Pour chaque point de la séquence, on calcule la p-value associée au test
  for (i in 1:(n-1)) {
    pvalues <- c(pvalues, test_func(x[1:i], x[(i+1):n])$p.value)
  }

  return(which(pvalues < threshold) + 1)
}

######################################################################
######################################################################

#' Seconde approche naïve
#'
#' @description Connaissant le nombre de ruptures, cette fonction segemente astucieusement pour la recherche de points de rupture avec, au choix le test de Wilcoxon ou de Kologorov-Smirnov.
#'
#' @details A chaque itération, on cherche l'indice associé à la plus petite p-value obtenue dans les sous-séquences données à l'itération précédente.
#' Attention, on suppose l'existence d'un point de rupture au minimum.
#'
#' @param x vecteur (série temporelle) à valeurs réelles
#' @param K supposé nombre de ruptures
#' @param method statistique à utiliser parmi "wilcox" et "KS"
#' @return \code{index_break} les indices associées aux ruptures
naive_test2 <- function(x, K = 1, method = c("wilcox", "KS")) {

  n <- length(x)

  if (method=="wilcox") { test_func <- wilcoxtest }
  if (method=="KS") { test_func <- ks.test }

  index_break <- NULL
  pvalues <- NULL

  # Recherche du 1er point de rupture le plus vraisemblable
  for (i in 1:(n-1)) {
    pvalues <- c(pvalues, test_func(x[1:i], x[(i+1):n])$p.value)
  }
  index_break <- c(index_break, which.min(pvalues) + 1)

  k <- 1
  # Recherche des autres points de rupture
  while (k < K) {
    pvalues <- NULL

    # Recherche à l'intérieur de la première composante à l'instant k
    for (i in 1:(index_break[1]-2)) {
      pvalues <- c(pvalues, test_func(x[1:i], x[(i+1):(index_break[1]-1)])$p.value)
    }

    # Recherche à l'intérieur des autres composantes
    if (k > 1) {
      for (j in 1:(k-1)) {
        for (i in (index_break[j]+1):(index_break[j+1]-2)) {
          pvalues <- c(pvalues, test_func(x[index_break[j]:i], x[(i+1):index_break[j+1]-1])$p.value)
        }
      }
    }

    # Recherche à l'intérieur de la dernière composante à l'instant k
    for (i in (index_break[k]+1):(n-1)) {
      pvalues <- c(pvalues, test_func(x[index_break[k]:i], x[(i+1):n])$p.value)
    }

    index_break <- c(index_break, which.min(pvalues) + 1)
    k <- k + 1
  }

  return(index_break)
}

######################################################################
######################################################################

#' Détection d'une rupture avec CUSUM
#'
#' @description Elle permet de détecter un point de rupture avec CUSUM sous l'hypothèse que l'on connaît la moyenne avant la rupture.
#'
#' @details Cette fonction va calculer des sommes partielles et retourner l'indice de rupture tel que la statistique dépasse un certain seuil en premier.
#'
#' @param x vecteur (série temporelle) à valeurs réelles
#' @param threshold seuil pour décider d'un changement de tendance
#' @param mu0 moyenne avant la rupture
#' @return \code{cp} le supposé indice du point de rupture \code{Tn} les statistiques de CUSUM
#' @references Romano, Gaetano and Eckley, Idris and Fearnhead, Paul and Rigaill, Guillem (2021) Fast Online Changepoint Detection via Functional Pruning CUSUM statistics, arXiv, 10.48550/ARXIV.2110.08205
#' @seealso https://github.com/gtromano/FOCuS
CUSUM <- function(x, threshold, mu0) {

  # Calcul des statistiques de CUSUM récursivement
  Tn <- abs(cumsum(x - mu0)) * (1/sqrt(1:length(x)))

  # Décision de rupture ?
  # Le premier qui a une statistique supérieure au seuil
  cp <- which(Tn >= threshold)[1]
  cp <- ifelse(is.na(cp), -1, cp)

  return(list(cp = cp, maxs = Tn))
}

######################################################################
######################################################################

#' Détection d'une rupture avec Page-CUSUM
#'
#' @description Elle permet de détecter un point de rupture avec Page-CUSUM sous l'hypothèse que l'on connaît la moyenne avant la rupture.
#'
#' @details Cette fonction va calculer des sommes partielles et retourner l'indice de rupture tel que la statistique dépasse un certain seuil en premier.
#'
#' @param x vecteur (série temporelle) à valeurs réelles
#' @param threshold seuil pour décider d'un changement de tendance
#' @param mu0 moyenne avant la rupture
#' @return \code{cp} le supposé indice du point de rupture et \code{Tn} les statistiques de Page-CUSUM
#' @references Romano, Gaetano and Eckley, Idris and Fearnhead, Paul and Rigaill, Guillem (2021) Fast Online Changepoint Detection via Functional Pruning CUSUM statistics, arXiv, 10.48550/ARXIV.2110.08205
#' @seealso https://github.com/gtromano/FOCuS
pageCUSUM <- function (x, threshold, mu0) {

  # Calcul des statistiques de Page-CUSUM récursivement
  Tn <- sapply(1:length(x), function (n) {
    Q <- abs(cumsum(x[n:1] - mu0)) * (1/sqrt(1:n))
    return(max(Q))
  })

  # Décision de rupture ?
  # Le premier qui a une statistique supérieure au seuil
  cp <- which(Tn >= threshold)[1]
  cp <- ifelse(is.na(cp), -1, cp)

  return(list(cp = cp, maxs = Tn))
}

######################################################################
######################################################################

#' Détection de ruptures multiples avec CUSUM ou Page-CUSUM
#'
#' @description Elle permet de détecter plusieurs points de rupture avec, au choix, CUSUM ou Page-CUSUM.
#'
#' @detials Lorsqu'un point de rupture est trouvé, la fonction va chercher sur les tronçons suivants d'autres potentiels points de rupture.
#' L'astuce utilisée est de remplacer les moyennes à priori (connues) par une donnée du vecteur \code{x}.
#' Attention, on suppose l'existence d'un point de rupture au minimum.
#'
#' @param x vecteur (série temporelle) à valeurs réelles
#' @param threshold seuil pour décider d'un changement de tendance
#' @param method statistique à utiliser parmi "CUSUM" et "Page-CUSUM"
#' @return \code{index_break} les indices associées aux ruptures
CUSUMs_test <- function(x, threshold, method = c("CUSUM", "Page-CUSUM")) {

  n <- length(x)
  index_break <- NULL

  if (method=="CUSUM") { CUSUM_func <- CUSUM }
  if (method=="Page-CUSUM") { CUSUM_func <- pageCUSUM }

  # Recherche du 1er point de rupture
  last_break <- CUSUM_func(x, threshold, x[1])$cp
  index_break <- c(index_break, last_break)

  # Recherche des potentiels autres points de rupture
  while(last_break<n) {
    temp_break <- CUSUM_func(x[(last_break+1):n], threshold, x[last_break])$cp
    if (temp_break==-1) { break } # Si plus aucun point de rupture

    last_break <- temp_break + last_break
    index_break <- c(index_break, last_break)
  }

  return(index_break)
}

######################################################################
######################################################################

#' Approche séquentielle de Page
#'
#' @description Elle permet de détecter un point de rupture séquentiellement sous l'hypothèse que l'on connaît la moyenne après la rupture.
#'
#' @details Cette fonction va calculer des rapports de log-vraisemblances et retourner l'indice de rupture tel que la statistique finale est non nulle (on suppose d'ailleurs que la moyenne avant la rupture est nulle).
#'
#' @param x vecteur (série temporelle) à valeurs réelles
#' @param threshold seuil pour décider d'un changement de tendance
#' @param mu1 moyenne après la rupture
#' @return \code{cp} le supposé indice du point de rupture et \code{Qn} les statistiques de la méthode séquentielle de Page (somme des rapports des vraisemblances logarithmiques)
pageSeq <- function(x, threshold, mu1) {

  # Calcul des statistiques de Page récursivement
  Qn <- cumsum(mu1 * (x - (mu1/2)))

  # Placer des 0 si la statistique est négative
  n <- length(x)
  for (i in 1:n) {
    if (Qn[i]<0) { Qn[i:n] <- Qn[i:n] - Qn[i]}
  }

  # Décision de rupture ?
  # Le point qui maximise la statistique
  cp <- which(Qn > threshold)[1]
  cp <- ifelse(is.na(cp), -1, cp)

  return(list(cp = cp, maxs = Qn))
}

######################################################################
######################################################################

#' Itération de l'algorithme FOCuS0
#'
#' @description Q_iter permet de mettre à jour Q avec un nouveau triplet (tau, s, l).
#'
#' @details Cette fonction va effectuer une mise à jour de Q par rapport à la nouvelle donnée pour l'algorithme de FOCuS0.
#'
#' @param Q matrice de triplets à valeurs réelles
#' @param s somme cumulative
#' @param n instant de la nouvelle observation
#' @return \code{Q} matrice après itération
#' @references Romano, Gaetano and Eckley, Idris and Fearnhead, Paul and Rigaill, Guillem (2021) Fast Online Changepoint Detection via Functional Pruning CUSUM statistics, arXiv, 10.48550/ARXIV.2110.08205
#' @seealso https://github.com/gtromano/FOCuS
Q_iter <- function(Q, s, n) {

  tau <- n
  l <- Inf
  k <- nrow(Q)
  i <- k

  # boucle while qui décrémente la valeur de i jusqu'à ce que la condition
  # (2 * (s - Q[i, 2]) - (tau - Q[i, 1]) * Q[i, 3] <= 0) soit fausse ou jusqu'à ce que i soit inférieur à 1
  while ((2 * (s - Q[i, 2]) - (tau - Q[i, 1]) * Q[i, 3] <= 0) && (k > 1) && (i >= 2)) {
    i <- i - 1
  }

  l <- max(0, 2 * (s - Q[i, 2]) / (tau - Q[i, 1]))

  # Si i < k, on supprime les lignes après les i premières lignes de Q
  if (i < k) { Q <- Q[1:i, ] }

  # Une nouvelle ligne est ajoutée à Q avec les valeurs tau, s et l.
  Q <- rbind(Q, c(tau, s, l))

  return(Q)
}

######################################################################
######################################################################

#' Détection d'une rupture avec FOCuS (R)
#'
#' @description FOCuS0 permet de détecter un point de rupture avec ou sans la connaissance de la moyenne après la rupture.
#'
#' @details Cette fonction va effectuer la première étape de l'algorithme FOCuS0, à savoir calculer Qn en fonction de différents mu1.
#' Puis (seconde étape) chercher la maximum sur tous les mu de Qn et enfin vérifier si la statistique est plus grand qu'un certain seuil.
#'
#' @param x vecteur (série temporelle) à valeurs réelles
#' @param threshold seuil pour décider d'un changement de tendance
#' @return \code{cp} le supposé indice du point de rupture et \code{Q} la matrice des itérations des triplets
#' @references Romano, Gaetano and Eckley, Idris and Fearnhead, Paul and Rigaill, Guillem (2021) Fast Online Changepoint Detection via Functional Pruning CUSUM statistics, arXiv, 10.48550/ARXIV.2110.08205
#' @seealso https://github.com/gtromano/FOCuS
FOCuS_R <- function(x, threshold) {

  n <- length(x)
  # Initialisation de Q (valant 0 nous n = 0)
  Q <- matrix(c(0, 0, 0), nrow = 1, ncol = 3)

  # Calculer la somme cumulée des éléments de x
  S <- cumsum(x)

  cp <- -1
  stop_F <- FALSE
  for (i in 1:n) {
    # Exécution de Q_iter avec les entrées Q et le i-ème élément de S
    Q <- Q_iter(Q, S[i], i)

    # On arrête l'algorithme lorsque la statistique est plus grande que threshold
    # pour l'un des quadratiques
    for (i_q in 1:nrow(Q)) {
      if ((S[i] - Q[i_q, 2])^2 >= 2 * threshold * (i - Q[i_q, 1])) {
        cp <- Q[i_q, 1]
        stop_F <- TRUE
        break
      }
    }

    if (stop_F) { break }
  }

  return (list(cp = cp, maxs = Q))
}
