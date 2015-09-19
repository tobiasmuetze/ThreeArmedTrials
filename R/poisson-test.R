#' @title Retention of effect hypothesis for poisson distributed endpoints
#' @description The function \code{test_RET_poisson} performs the Wald-type test
#' for the retention of effect hypothesis with poisson distributed endpoints
#' @details The hypothesis
#' \eqn{(\lambda_P - \lambda_E)/(\lambda_P - \lambda_R) \le \Delta}
#' is tested against the alternative
#' \eqn{(\lambda_P - \lambda_E)/(\lambda_P - \lambda_R) > \Delta}.
#' Here, \eqn{\lambda_E}, \eqn{\lambda_R}, \eqn{\lambda_P}
#' are the rates of the experimental treatment (\code{rateExp}), the reference
#' treatment (\code{rateRef}), and the placebo group (\code{ratePla}), respectively.
#' @param xExp A (non-empty) numeric vector of data values coming from the experimental treatment group.
#' @param xRef A (non-empty) numeric vector of data values coming from the reference treatment group.
#' @param xPla A (non-empty) numeric vector of data values coming from the placebo group.
#' @param Delta A numeric value specifying the margin.
#' @return A list with class "htest" containing the following components:
#' \item{statistic}{The value of the Wald-type test statistic.}
#' \item{p.value}{The p-value for the Wald-type test.}
#' \item{method}{A character string indicating what type of Wald-type-test was performed.}
#' \item{estimate}{The estimated rates for each of the group as well as the maximum-likelihood estimator for the shape parameter.}
#' \item{sample.size}{The total number of data points used for the Wald-type test.}
#' @examples
#' # Poisson distributed endpoints
#' # Test for non-inferiority test. lambda_P = 8, lambda_R = 4, lambda_E = 5
#' # Delta = (lambda_P-lambda_E)/(lambda_P-lambda_R)
#' xExp <- rpois(60, lambda = 5)
#' xRef <- rpois(40, lambda = 4)
#' xPla <- rpois(40, lambda = 8)
#' Delta <- (8-5)/(8-4)
#' test_RET_poisson(xExp, xRef, xPla, Delta)
#'
#' # Test for superiority test. lambda_P = 8, lambda_R = 5, lambda_E = 4
#' # Delta = (lambda_P-lambda_E)/(lambda_P-lambda_R)
#' xExp <-rpois(60, lambda = 5)
#' xRef <-rpois(40, lambda = 4)
#' xPla <-rpois(40, lambda = 8)
#' Delta <- (8-5)/(8-4)
#' test_RET_poisson(xExp, xRef, xPla, Delta)
#' @references Mielke, M. and Munk, A. (2010).
#' \emph{The assessment and planning of non-inferiority trials for retention of effect hypotheses - towards a general approach.}
#' arXiv:0912.4169v1
#' @export
#' @keywords test Poisson
test_RET_poisson <- function(xExp, xRef, xPla, Delta){

  # To-Do: Alternative 'smaller'

  if (!is.numeric(c(xExp, xRef, xPla))) {
    stop('Data must be numeric.')
  }

  if (!is.numeric(Delta) || (Delta < 0)) {
    stop('Delta must be postive.')
  }

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

  # Rate estimators
  rateExp <- mean(xExp)
  rateRef <- mean(xRef)
  ratePla <- mean(xPla)

  # variance for test statistic
  varWaldTest <- rateExp / wExp + Delta^2 * rateRef / wRef + (Delta - 1)^2 * ratePla / wPla

  # Test statistic and p-value
  teststat <- sqrt(n) * (  rateExp - Delta * rateRef + (Delta-1) * ratePla ) / sqrt(varWaldTest)
  pvalue <- pnorm(teststat)


  # Output
  estimate <- c(rateExp, rateRef, ratePla)
  names(estimate) <- c('Rate Exp', 'Rate Ref', 'Rate Pla')
  names(teststat) <- c('T')
  methodText <- 'Wald-type test with maximum-likelihood variance estimator'
  structure(list(statistic = teststat,
                 p.value = pvalue,
                 method = methodText,
                 data.name = data.name,
                 estimate = estimate,
                 sample.size = n),
            class = "htest")
}






#' @title Retention of effect hypothesis for poisson distributed endpoints
#' @description The function \code{simulate_test_RET_poisson} performs a
#' simulation study of the rejection rate of the Wald-type test
#' for the retention of effect hypothesis with poisson distributed endpoints
#' @details The hypothesis
#' \eqn{(\lambda_P - \lambda_E)/(\lambda_P - \lambda_R) \le \Delta}
#' is tested against the alternative
#' \eqn{(\lambda_P - \lambda_E)/(\lambda_P - \lambda_R) > \Delta}.
#' Here, \eqn{\lambda_E}, \eqn{\lambda_R}, \eqn{\lambda_P}
#' are the rates of the experimental treatment (\code{rateExp}), the reference
#' treatment (\code{rateRef}), and the placebo group (\code{ratePla}), respectively.
#' @param rateExp A numeric values defining the rate in the experimental treatment group.
#' @param rateRef A numeric values defining the rate in the reference treatment group.
#' @param ratePla A numeric values defining the rate in the placebo treatment group.
#' @param nExp A numeric values defining the sample size in the  experimental treatment group.
#' @param nRef A numeric values defining the rate in the reference treatment group.
#' @param nPla A numeric values defining the rate in the placebo treatment group.
#' @param Delta A positive numeric value specifying the margin.
#' @param nSimulations Number of Monte-Carlo samples used for simulation study.
#' @param sig.level A numeric value specifying the significance level for the Wald-tye test.
#' @return A list with class "htest" containing the following components:
#' \item{statistic}{The value of the Wald-type test statistic.}
#' \item{p.value}{The p-value for the Wald-type test.}
#' \item{method}{A character string indicating what type of Wald-type-test was performed.}
#' \item{estimate}{The estimated rates for each of the group as well as the maximum-likelihood estimator for the shape parameter.}
#' \item{sample.size}{The total number of data points used for the Wald-type test.}
#' @examples
#' # Poisson distributed endpoints
#' # Test for non-inferiority test. lambda_P = 8, lambda_R = 4, lambda_E = 5
#' # Delta = (lambda_P-lambda_E)/(lambda_P-lambda_R)
#' Delta <- (8-5)/(8-4)
#' @references Mielke, M. and Munk, A. (2010).
#' \emph{The assessment and planning of non-inferiority trials for retention of effect hypotheses - towards a general approach.}
#' arXiv:0912.4169v1
#' @export
#' @keywords test Poisson
simulate_test_RET_poisson <- function(rateExp, rateRef, ratePla, nExp, nRef, nPla, nSimulations, Delta, sig.level){

  if (!is.numeric(Delta) || (Delta < 0)) {
    stop('Delta must be postive.')
  }
  if (!is.numeric(sig.level) || (sig.level <= 0) || (sig.level >= 1)) {
    stop('Significance level must be between 0 and 1.')
  }

  n <- nExp + nRef + nPla

  wExp <- nExp / n
  wRef <- nRef / n
  wPla <- nPla / n

  xExp <- matrix( rpois(nExp * nSimulations, lambda = rateExp), ncol = nExp, nrow = nSimulations)
  xRef <- matrix( rpois(nRef * nSimulations, lambda = rateRef), ncol = nRef, nrow = nSimulations)
  xPla <- matrix( rpois(nPla * nSimulations, lambda = ratePla), ncol = nPla, nrow = nSimulations)



  meanExp <- rowMeans(xExp)
  meanRef <- rowMeans(xRef)
  meanPla <- rowMeans(xPla)

  varWaldTest <- meanExp / wExp + Delta^2 * meanRef / wRef + (Delta - 1)^2 * meanPla / wPla

  # Test statistic and p-value
  teststat <- sqrt(n) * (  meanExp - Delta * meanRef + (Delta-1) * meanPla ) / sqrt(varWaldTest)
  pvalue <- pnorm(teststat)


  rejection_rate <- mean(pvalue <= sig.level)

  data.frame(rateExp, rateRef, ratePla, nExp, nRef, nPla, nSimulations, Delta, sig.level, rejection_rate)

}
