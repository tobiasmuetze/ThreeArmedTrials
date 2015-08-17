#' @title Studentized permutation test analyzing for three-arm clinical trials.
#' @description Wald-type test for superiority/non-inferiority of the experimental treatment versus reference treatment with respect to placebo.
#' @details The hypothesis \eqn{(\lambda_P - \lambda_E)/(\lambda_P - \lambda_R) \le \Delta} is tested against the alternative \eqn{(\lambda_P - \lambda_E)/(\lambda_P - \lambda_R) > \Delta}.
#' \eqn{\lambda_E}, \eqn{\lambda_R}, \eqn{\lambda_P} are the expected responses of the experimental treatment (\code{rateExp}), the reference treatment (\code{rateRef}), and the placebo group (\code{ratePla}), respectively.
#' The margin \emph{Delta}, i.e. \eqn{\Delta} in the formulas above, is between 0 and 1 for testing non-inferiority and larger than 1 for testing superiority.
#' @param xExp A (non-empty) numeric vector of data values coming from the experimental treatment group.
#' @param xRef A (non-empty) numeric vector of data values coming from the reference treatment group.
#' @param xPla A (non-empty) numeric vector of data values coming from the placebo group.
#' @param Delta A numeric value specifying the non-inferiority or superiority margin. Is between 0 and 1 in case of non-inferiority and larger than 1 in case of superiority.
#' @param nPermutations Number of permutations.
#' @return A list with class "htest" containing the following components:
#' \item{statistic}{The value of the Wald-type test statistic.}
#' \item{p.value}{The p-value for the Wald-type test.}
#' \item{method}{A character string indicating what type of Wald-type-test was performed.}
#' \item{estimate}{The estimated rates for each of the group as well as the maximum-likelihood estimator for the shape parameter.}
#' \item{parameters}{Margin Delta.}
#' \item{sample.size}{The total number of data points used for the Wald-type test.}
#' \item{nPermutations}{Permutation distribution of the test statistic.}
#' @usage permutation.test(xExp, xRef, xPla, Delta, nPermutations)
#' @examples permutation.test(xExp = rnorm(15, mean=4),
#' xPla = rnorm(15, mean=5), xRef = rnorm(15, mean=3),
#' Delta = 0.5, nPermutations = 5000)
#' @export
#' @keywords test permutation
permutation.test <- function(xExp, xRef, xPla, Delta, nPermutations){

  if(Delta<=0){stop('Margin \'Delta\' must be postive.')}
  if( (nPermutations<=0) || !(nPermutations==round(nPermutations)) ){stop('Number of permutations \'nPermutations\' must be postive whole number.')}

  if(!is.numeric(c(xExp, xRef, xPla))){stop('Input data must be numeric.')}

  data.name <- paste(c(deparse(substitute(xExp)), ', ', deparse(substitute(xRef)), ', and ', deparse(substitute(xPla)), collapse=''))

  xExp <- xExp[!is.na(xExp)]
  xRef <- xRef[!is.na(xRef)]
  xPla <- xPla[!is.na(xPla)]

  # Sample size allocation
  nExp <- length(xExp)
  nRef <- length(xRef)
  nPla <- length(xPla)
  n <- nExp + nRef + nPla
  wExp <- nExp / n
  wRef <- nRef / n
  wPla <- nPla / n

  # Test statistic
  sigma2.Teststat <- var(xExp) / wExp + Delta^2 * var(xRef) / wRef + (1-Delta)^2 * var(xPla) / wPla
  Teststat <- sqrt(n) * (  mean(xExp) - Delta * mean(xRef) + (Delta-1) * mean(xPla) ) / sqrt(sigma2.Teststat)


  # Permutation test
  x <- c(xExp, xRef, xPla)
  n <- length(x)
  xPerm <- t(sapply(X = 1:nPermutations, function(i){sample(x, replace = F)} ))

  xExpPerm <- xPerm[,1:nExp]
  xRefPerm <- xPerm[,(nExp+1):(nExp+nRef)]
  xPlaPerm <- xPerm[,(nExp+nRef+1):n]

  xExpPermMean <- rowMeans(xExpPerm)
  xRefPermMean <- rowMeans(xRefPerm)
  xPlaPermMean <- rowMeans(xPlaPerm)

  sigma2ExpEst <- ( rowSums(xExpPerm^2) - nExp * xExpPermMean^2 ) / (nExp - 1)
  sigma2RefEst <- ( rowSums(xRefPerm^2) - nRef * xRefPermMean^2 ) / (nRef - 1)
  sigma2PlaEst <- ( rowSums(xPlaPerm^2) - nPla * xPlaPermMean^2 ) / (nPla - 1)

  sigma2.TeststatPerm <- sigma2ExpEst / wExp + Delta^2 * sigma2RefEst / wRef + (1-Delta)^2 * sigma2PlaEst / wPla
  TeststatPerm <- sqrt(n) * (  xExpPermMean - Delta * xRefPermMean + (Delta-1) * xPlaPermMean ) / sqrt(sigma2.TeststatPerm)

  pValuePerm <- sum( Teststat >= TeststatPerm ) / (nPermutations+1) + 1/ (nPermutations+1)

  # Output
  estimate <- c(mean(xExp), mean(xRef), mean(xPla))
  names(estimate) <- c('Mean Exp', 'Mean Ref', 'Mean Pla')
  names(Teststat) <- c('T')
  names(Delta) <- 'Delta'
  methodText <- paste0('Studentized permutation test with ', nPermutations, ' permutations')
  out <- structure(list(statistic = Teststat,
                 p.value = pValuePerm,
                 parameters = Delta,
                 method = methodText,
                 data.name = data.name,
                 estimate = estimate,
                 sample.size = n,
                 permDist = TeststatPerm),
            class = "htest")

return(out)
}





#' @title Simulation of reject rate for studentized permutation test within three-arm trials.
#' @description Wald-type test for superiority/non-inferiority of the experimental treatment versus reference treatment with respect to placebo.
#' @details The hypothesis \eqn{(\lambda_P - \lambda_E)/(\lambda_P - \lambda_R) \le \Delta} is tested against the alternative \eqn{(\lambda_P - \lambda_E)/(\lambda_P - \lambda_R) > \Delta}.
#' \eqn{\lambda_E}, \eqn{\lambda_R}, \eqn{\lambda_P} are the expected responses of the experimental treatment (\code{rateExp}), the reference treatment (\code{rateRef}), and the placebo group (\code{ratePla}), respectively.
#' The margin \emph{Delta}, i.e. \eqn{\Delta} in the formulas above, is between 0 and 1 for testing non-inferiority and larger than 1 for testing superiority.
#' @param parameterExp A (non-empty) numeric vector of data values coming from the experimental treatment group.
#' @param distExp A (non-empty) numeric vector of data values coming from the reference treatment group.
#' @param parameterRef A (non-empty) numeric vector of data values coming from the placebo group.
#' @param distRef A numeric value specifying the non-inferiority or superiority margin. Is between 0 and 1 in case of non-inferiority and larger than 1 in case of superiority.
#' @param parameterPla A (non-empty) numeric vector of data values coming from the placebo group.
#' @param distPla A numeric value specifying the non-inferiority or superiority margin. Is between 0 and 1 in case of non-inferiority and larger than 1 in case of superiority.
#' @param Delta A numeric value specifying the non-inferiority or superiority margin. Is between 0 and 1 in case of non-inferiority and larger than 1 in case of superiority.
#' @param nPermutations Number of permutations for permutation test.
#' @param nSimulations Number of Monte-Carlo samples used for simulation study.
#' @param sig.level A numeric value specifying the significance level for the permutation test.
#' @usage simulation.permutation.test(parameterExp, distExp, parameterRef, distRef,
#' parameterPla, distPla, nSimulations, nPermutations, Delta, sig.level)
#' @examples # Simulate significance level when data is negative binomially distributed
#' simulation.permutation.test(parameterExp=list(n=30,mu=7.56,size=1),distExp='rnbinom',
#' parameterRef=list(n=30,mu=5.1, size=1), distRef='rnbinom',
#' parameterPla=list(n=30,mu=17.4, size=1), distPla='rnbinom',
#' nSimulations=500, nPermutations=1000, Delta=0.8, sig.level=0.025)
#' @export
#' @keywords test permutation simulation
simulation.permutation.test <- function(parameterExp, distExp, parameterRef, distRef, parameterPla, distPla, nSimulations, nPermutations, Delta, sig.level){

  reject <- 0

  # get total sample size
  xExp <- do.call(what = distExp, args = parameterExp)
  xRef <- do.call(what = distRef, args = parameterRef)
  xPla <- do.call(what = distPla, args = parameterPla)
  nExp <- length(xExp)
  nRef <- length(xRef)
  nPla <- length(xPla)
  n <- nExp + nRef + nPla
  wExp <- nExp / n
  wRef <- nRef / n
  wPla <- nPla / n

  # Initialize permutations
  permMat <- t(sapply(X = 1:nPermutations, function(i){sample(1:n, replace = F)} ))

  for(i in 1:nSimulations){
    xExp <- do.call(what = distExp, args = parameterExp)
    xRef <- do.call(what = distRef, args = parameterRef)
    xPla <- do.call(what = distPla, args = parameterPla)

    # Test statistic for current data set
    sigma2.Teststat <- var(xExp) / wExp + Delta^2 * var(xRef) / wRef + (1-Delta)^2 * var(xPla) / wPla
    Teststat <- sqrt(n) * (  mean(xExp) - Delta * mean(xRef) + (Delta-1) * mean(xPla) ) / sqrt(sigma2.Teststat)

    # Permutation data
    xPerm <- matrix(c(xExp, xRef, xPla)[permMat], ncol=n)

    xExpPerm <- xPerm[,1:nExp]
    xRefPerm <- xPerm[,(nExp+1):(nExp+nRef)]
    xPlaPerm <- xPerm[,(nExp+nRef+1):n]

    xExpPermMean <- rowMeans(xExpPerm)
    xRefPermMean <- rowMeans(xRefPerm)
    xPlaPermMean <- rowMeans(xPlaPerm)

    sigma2ExpEst <- ( rowSums(xExpPerm^2) - nExp * xExpPermMean^2 ) / (nExp - 1)
    sigma2RefEst <- ( rowSums(xRefPerm^2) - nRef * xRefPermMean^2 ) / (nRef - 1)
    sigma2PlaEst <- ( rowSums(xPlaPerm^2) - nPla * xPlaPermMean^2 ) / (nPla - 1)

    sigma2.TeststatPerm <- sigma2ExpEst / wExp + Delta^2 * sigma2RefEst / wRef + (1-Delta)^2 * sigma2PlaEst / wPla
    TeststatPerm <- sqrt(n) * (  xExpPermMean - Delta * xRefPermMean + (Delta-1) * xPlaPermMean ) / sqrt(sigma2.TeststatPerm)

    pValue <- sum( Teststat >= TeststatPerm ) / (nPermutations+1) + 1/ (nPermutations+1)

    if(pValue <= sig.level ){
      reject <- reject + 1
    }
  }

  rejectRate <- reject / nSimulations
  names(parameterExp) <- paste0(names(parameterExp), "_exp")
  names(parameterRef) <- paste0(names(parameterRef), "_ref")
  names(parameterPla) <- paste0(names(parameterPla), "_pla")

  out <- c(parameterExp, parameterRef, parameterPla, list(rejectRate=rejectRate))
  return(as.data.frame(out))
}


