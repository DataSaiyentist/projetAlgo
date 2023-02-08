######################################################################
###                                                                ###
###  GPL-3 License                                                 ###
###  Copyright (c) 2023 @gtromano                                  ###
###                                                                ###
######################################################################

#' @title Fast Online Changepoint Detection via Functional Pruning CUSUM statistics (C++)
#'
#' @description FOCuS est un algorithme permettant de détecter les changements de moyenne en temps réel. Ceci est réalisé par une mise à jour récursive d'une équation du second degré par morceaux, dont le maximum est la statistique de test CUSUM pour un changement. FOCuS peut être appliqué à des situations où la moyenne avant changement est connue ou inconnue. En outre, FOCuS peut détecter des changements en présence de points aberrants.
#'
#' @param datasource une fonction qui génère des données (cf. \strong{Details})
#' @param thres seuil pour décider d'un changement de tendance
#' @param mu0 moyenne avant la rupture, si connue. Par défaut est \code{NA} et dans ce cas, la moyenne avant changement est estimée itérativement
#' @param grid vecteur de valeurs des amplitudes de changement pour activer l'approximation de la grille de FOCuS. La valeur par défaut est \code{NA}, de sorte que, par défaut, FOCuS fonctionne exactement (cf. \strong{Details})
#' @param K La valeur de la fonction de coût bi-poids. La valeur par défaut est \code{Inf}. FOCuS ignorera les points aberrants plus grands que \code{K}
#'
#' @details \code{datasource} requiert un objet de la classe \linkS4class{function}. Une fonction valide doit fournir un point des données à chaque appel, et ne doit prendre aucun argument supplémentaire. À chaque itération, FOCuS extrait une observation par un appel de fonction, traite la nouvelle observation et répète cette procédure jusqu'à trouver un point de changement. (cf. \strong{Examples}).
#'
#' @return un objet s3 de la classe FOCuSout où :
#' \describe{
#' \item{\code{$t}}{est l'indice de détection,}
#' \item{\code{$changepoint}}{est le point de rupture,}
#' \item{\code{$Q1}}{est le coût optimal sous forme de d'une équation du second degré par morceaux à la fin de la séquence au moment de la détection.}
#' }
#'
#' @references Romano, Gaetano and Eckley, Idris and Fearnhead, Paul and Rigaill, Guillem (2021) Fast Online Changepoint Detection via Functional Pruning CUSUM statistics, arXiv, 10.48550/ARXIV.2110.08205
#' @seealso https://github.com/gtromano/FOCuS (Propriétaires de l'algorithme)
#'
#' @examples
#' ##################
#' ###   online   ###
#' ##################
#'
#' set.seed(42)
#' databuffer <- c(rnorm(3e5, 1), rnorm(1e4, 0))
#'
#' f <- function() {
#'   out <- databuffer[i]     # simulating a pull from a buffer
#'   i <<- i + 1
#'   out
#' }
#'
#' i <- 1; FOCuS(f, 18)

setGeneric(
  "FOCuS",
  def = function(datasource, thres, ...) standardGeneric("FOCuS")
)

setMethod("FOCuS",
          signature(data = "function", thres = "numeric"),
          function (datasource, thres,  mu0 = NA, grid = NA, K = Inf)
          {
            # checks on the data generating function
            if(is.function(match.fun(datasource)))
              out_check <- datasource()
            else
              stop("Please provvide a data generating function.")
            if(!is.numeric(out_check)) stop("Data generating function provvided does not return a numeric output.")
            if(length(out_check) != 1) {
              warning("Length of the output from data generating function is greater than 1. Using only the first element.")
              f <- function() return(datasource()[1])
            } else {
              f <- match.fun(datasource)
            }

            # checks on the thres
            if( !is.numeric(thres) | thres < 0)
              stop("thres must be a positive numeric")

            # checks on the mu0
            if(!is.na(mu0))
              if(!is.numeric(mu0) | length(mu0) > 1)
                stop("mu0 must be a numeric value")

            # checks on the grid
            if(!is.na(grid))
              if(!is.numeric(grid))
                stop("mu0 must be a numeric value")

            # checks on the K
            if(!is.na(K))
              if(!is.numeric(K) | K <= 0)
                stop("K must be a positive numeric")


            out <- .FoCUS(datasource, thres, mu0, grid, K)
            out$changepoint <- out$t + out$Q1[[which.max(sapply(out$Q1, function(q) q$max))]]$a * 2
            class(out) <-  c("FOCuSout", class(out))
            class(out$Q1) <- c("PiecewiseQuadratic", class(out))

            if (!is.null(out$warning_message))
              warning(out$warning_message)

            return(out)
          }
)
