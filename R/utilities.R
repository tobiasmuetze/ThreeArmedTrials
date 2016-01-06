#' @title rt_nonstandard
#' @description random numbers of the non-standard t-distribution
#' @param n number of observations.
#' @param df degrees of freedom (> 0, maybe non-integer). df = Inf is allowed.
#' @param mu Location parameter.
#' @param sigma Shape parameter.
rt_nonstandard <- function(n, df, mu, sigma) {
  rt(n = n, df = df) * sigma + mu
}

#' @title loglikelihood_binary
#' @description log likelihood of Bernoulli function
#' @param p numeric vector of probabilities with length 3
#' @param xExp numeric vector of probabilities with length 3
#' @param xRef numeric vector of probabilities with length 3
#' @param xPla numeric vector of probabilities with length 3
loglikelihood_binary <- function(p, xExp, xRef, xPla) {
  pExp <- p[1]
  pRef <- p[2]
  pPla <- p[3]

  nExp <- length(xExp)
  nRef <- length(xRef)
  nPla <- length(xPla)

  sumExp <- sum(xExp)
  sumRef <- sum(xRef)
  sumPla <- sum(xPla)

  log_l <- sumExp * log(pExp) + (nExp - sumExp) * log(1 - pExp) +
    sumRef * log(pRef) + (nRef - sumRef) * log(1 - pRef) +
    sumPla * log(pPla) + (nPla - sumPla) * log(1 - pPla)

  -log_l
}

#' @title is.naturalnumber
#' @description check if input is natural number
#' @param x numeric number to be checked
#' @param tol maximum accepted tolerance when checking if natural
is.naturalnumber <- function(x, tol = .Machine$double.eps^0.5) {
  ((abs(x - round(x)) < tol) & (x > 0))
}
