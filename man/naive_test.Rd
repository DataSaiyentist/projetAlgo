% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/naive_approach.R
\name{naive_test}
\alias{naive_test}
\title{Première approche naïve}
\usage{
naive_test(x, threshold = 0.05, method = c("wilcox", "KS"))
}
\arguments{
\item{x}{vecteur (série temporelle) à valeurs réelles}

\item{threshold}{seuil pour décider d'un changement de tendance}

\item{method}{statistique à utiliser parmi "wilcox" et "KS"}
}
\value{
\code{index_break} les indices associées aux ruptures
}
\description{
Cette fonction teste sur chaque point (votre choix) le test de Wilcoxon ou de Kologorov-Smirnov.
}
\details{
Elle retourne les indices des points pour lesquelles la p-value est plus petite q'un certain seuil.
Attention, le test de Kolmogorov-Smirnov n'a pas été implémenté à la main (\code{ks.test} de la librairie \code{"stats"}).
}
