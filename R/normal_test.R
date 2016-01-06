#' @title Statistical test for three-armed clinical trials with normally distributed endpoints.
#' @description Wald-type test for superiority/non-inferiority of the experimental treatment versus reference treatment with respect to placebo.
#' @details
#' The hypothesis
#' \eqn{(\lambda_P - \lambda_E)/(\lambda_P - \lambda_R) \le \Delta}
#' is tested against the alternative
#' \eqn{(\lambda_P - \lambda_E)/(\lambda_P - \lambda_R) > \Delta}.
#' \eqn{\lambda_E}, \eqn{\lambda_R}, \eqn{\lambda_P}
#' are the expected values of the normally distribtued endpoints in the
#' experimental treatment (\code{rateExp}), the reference treatment
#' (\code{rateRef}), and the placebo group (\code{ratePla}), respectively.
#' The variances \eqn{\sigma^2_{E}}, \eqn{\sigma^2_{R}},
#' and \eqn{\sigma^2_{P}} of the normally distributed
#' endpoint are by default heterogeneous.
#' The margin \emph{Delta}, i.e. \eqn{\Delta} in the
#' formulas above, is between 0 and 1 for testing
#' non-inferiority and larger than 1 for testing superiority.
#' @param xExp A (non-empty) numeric vector of data values coming from the experimental treatment group.
#' @param xRef A (non-empty) numeric vector of data values coming from the reference treatment group.
#' @param xPla A (non-empty) numeric vector of data values coming from the placebo group.
#' @param Delta A numeric variable specifying the non-inferiority or superiority margin.
#' Is between 0 and 1 in case of non-inferiority and larger than 1 in case of superiority.
#' @param var.equal A logical variable indicating whether to treat the three variances as being equal.
#' If TRUE then the pooled variance is used to estimate the variance
#' otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.
#' @return A list with class "htest" containing the following components:
#' \item{statistic}{The value of the t-statistic.}
#' \item{p.value}{The p-value for the Wald-type test.}
#' \item{var.equal}{A logical variable indicating what the three variances were treated as equal.}
#' \item{estimate}{The estimated means and variances for each of the groups.}
#' \item{parameter}{The degrees of freedom for the t-statistic.}
#' \item{sample.size}{The total number of data points used for the Wald-type test.}
#' \item{bad.obs}{The number of missing (NA), undefined (NaN) and/or infinite (Inf, -Inf) values that were removed from the input data objects prior to performing the hypothesis test.}
#' @examples
#' # Normally distributed endpoints with equal variances
#' # Test for non-inferiority test
#' xExp <-rnorm(60, mean=5, sd=1)
#' xRef <-rnorm(40, mean=4, sd=1)
#' xPla <-rnorm(40, mean=8, sd=1)
#' Delta <- (8-5)/(8-4)
#' test_RET_normal(xExp, xRef, xPla, Delta, var.equal = TRUE)
#'
#' # Normally distributed endpoints with heteroscedastic variances
#' # Test for non-inferiority test
#' xExp <-rnorm(60, mean=5, sd=1)
#' xRef <-rnorm(40, mean=4, sd=2)
#' xPla <-rnorm(40, mean=8, sd=3)
#' Delta <- (8-5)/(8-4)
#' test_RET_normal(xExp, xRef, xPla, Delta, var.equal = FALSE)
#' @references Hasler, M., Vonk, R., and Hothorn, L.A. (2008).
#' \emph{Assessing non-inferiority of a new treatment in a three-arm trial in the presence of heteroscedasticity.}
#' Statistics in Medicine 27, 490-503.
#'
#' Pigeot, I., Schaefer, J., Roehmel, J., and Hauschke, D. (2003).
#' \emph{Assessing non-inferiority of a new treatment in a three-arm clinical trial including a placebo.}
#' Statistics in Medicine 22, 883-899.
#' @seealso \code{\link{power.taNegbin.test}}
#' @export
#' @keywords test NormalDistribution
test_RET_normal <- function(xExp, xRef, xPla, Delta, var.equal = FALSE){

  data.name <- paste(c(deparse(substitute(xExp)), ', ',
                       deparse(substitute(xRef)), ', and ',
                       deparse(substitute(xPla)), collapse=''))

  # var.equal
  if (!is.logical(var.equal)) {
    stop("'var.equal' is not a logical variable.")
  }
  # Number of bad observations
  bad.obs <- sum(is.na(c(xExp, xRef, xPla)))
  # NA entries
  if ( any(is.na(c(xExp, xRef, xPla)))) {
    warning("'xExp', 'xRef', or 'xPla' contain NA entries which will be removed.")
    xExp <- xExp[!is.na(xExp)]
    xRef <- xRef[!is.na(xRef)]
    xPla <- xPla[!is.na(xPla)]
  }
  # Non-numeric entries
  if (any(!is.numeric(c(xExp, xRef, xPla)))) {
    stop("'xExp', 'xRef', must 'xPla' contain only numeric entries.")
  }
  # Delta
  if (!is.numeric(Delta)) {
    stop("Margin 'Delta' must be a numeric variable.")
  }
  if (Delta <= 0) {
    stop("Margin 'Delta' must be postive.")
  }

  # Sample size allocation
  nExp <- length(xExp)
  nRef <- length(xRef)
  nPla <- length(xPla)

  if (any(c(nExp, nRef, nPla) < 2)) {
    stop("'xExp', 'xRef', and 'xPla' must each have a length of at least 2.")
  }

  n <- nExp + nRef + nPla
  wExp <- nExp / n
  wRef <- nRef / n
  wPla <- nPla / n

  # Rate estimators
  meanExp <- mean(xExp)
  meanRef <- mean(xRef)
  meanPla <- mean(xPla)

  # Variance estimator
  if (var.equal) {
    varWaldTest <- (1/nExp + Delta^2 / nRef + (1-Delta)^2 / nPla) *
      ( (nExp-1)*var(xExp) + (nRef-1)*var(xRef) + (nPla-1)*var(xPla) ) / (n-3)
    degF <- n-3
  } else if (!var.equal) {
    varWaldTest <- var(xExp) / nExp +
      Delta^2 * var(xRef) / nRef +
      (1-Delta)^2 * var(xPla) / nPla
    if (varWaldTest == 0) {
      degF <- 1
    } else {
      degF <- varWaldTest^2 /
        ( var(xExp)^2 / (nExp^2 * (nExp-1)) +
            Delta^4 * var(xRef)^2 / (nRef^2 * (nRef-1)) +
            (1-Delta)^4 * var(xPla)^2 / (nPla^2 * (nPla-1))  )
    }
  }


  # Test statistic and p-value
  effect_estimator <- meanExp - Delta * meanRef + (Delta-1) * meanPla
  if ((effect_estimator == 0) && (varWaldTest == 0)) {
    Teststat <- 0
  } else {
    Teststat <- effect_estimator / sqrt(varWaldTest)
  }
  pvalue <- pt(Teststat, df = degF)

  if (all(c(xExp, xRef, xPla) == 0)) {
    pvalue <- 1
  }

  # Output
  estimate <- c(meanExp, meanRef, meanPla)
  names(estimate) <- c('Mean Exp', 'Mean Ref', 'Mean Pla')
  names(Teststat) <- c('T')
  if(var.equal){
    methodText <- 'Wald-type test under assumption of homogeneous variances'
  } else{
    methodText <- 'Wald-type test under assumption of heterogeneous variances'
  }

  structure(list(statistic = Teststat,
                 p.value = pvalue,
                 method = methodText,
                 data.name = data.name,
                 estimate = estimate,
                 parameter = degF,
                 sample.size = n,
                 bad.obs = bad.obs
                 ),
            class = "htest")
}
