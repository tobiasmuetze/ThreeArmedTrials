#' @title Power related calcuations for retention of effect hypothesis with Poisson distributed endpoints
#' @description Compute power, sample size, or significance level for Wald-type test
#' testing the retention of effect hypothesis with Poisson distributed endpoints.
#' @details If the individual group sample sizes, i.e. \code{n*allocation} are not natural number,
#' the parameters \emph{n} and \emph{allocation} will be re-calculated.
#' @param rateExp A numeric value specifying the rate of the experimental treatment group in the alternative hypothesis
#' @param rateRef A numeric value specifying the rate of the reference treatment group in the alternative hypothesis
#' @param ratePla A numeric value specifying the rate of the placebo treatment group in the alternative hypothesis
#' @param Delta A numeric value specifying the non-inferiority or superiority margin.
#' @param sig.level A numeric value specifying the significance level (type I error probability)
#' @param power A numeric value specifying the target power (1 - type II error probability)
#' @param n The total sample size. Needs to be at least 7.
#' @param allocation A (non-empty) vector specifying the sample size allocation (nExp/n, nRef/n, nPla/n)
#' @param type A character string determing how the variance for the Wald-type test statistic is estimated, \emph{unrestricted}
#' @return A list with class "power.htest" containing the following components:
#' \item{n}{The total sample size}
#' \item{power}{A numeric value specifying the target power}
#' \item{Delta}{A numeric value specifying the non-inferiority or superiority margin. }
#' \item{sig.level}{A character string specifying the significance level}
#' \item{type}{A character string indicating what type of Wald-type test will be performed}
#' \item{allocation}{A vector with the sample size allocation (nExp/n, nRef/n, nPla/n)}
#' \item{sig.level}{The significance level (Type I error probability)}
#' \item{nExp}{A numeric value specifying the number of sample in the experimental treatment group}
#' \item{nRef}{A numeric value specifying the number of sample in the reference treatment group}
#' \item{nPla}{A numeric value specifying the number of sample in the placebo treatment group}
#' @examples
#' # Initiate rates for defining margin in null hypothesis
#' rateExp <- 0.29
#' rateRef <- 0.2
#' ratePla <- 0.40
#'
#' # Margin in the hypothesis
#' Delta <- (ratePla - rateExp) / (ratePla - rateRef)
#'
#' # Define equal rates in the alternative
#' rateExp <- rateRef
#'
#' sig.level <- 0.025
#' power_level <- 0.8
#'
#' # Determine Sample Size
#' power_RET_poisson(rateExp = rateExp,
#'                   rateRef = rateRef,
#'                   ratePla = ratePla,
#'                   Delta = Delta,
#'                   sig.level = sig.level,
#'                   power = power_level,
#'                   type = 'unrestricted',
#'                   allocation = c(1/3, 1/3, 1/3))
#'
#' # Determine Power
#' power_RET_poisson(rateExp = rateExp,
#'                  rateRef = rateRef,
#'                  ratePla = ratePla,
#'                  Delta = Delta,
#'                  sig.level = sig.level,
#'                  n = 993,
#'                  type = 'unrestricted',
#'                  allocation = c(1/3, 1/3, 1/3))
#' @references Mielke, M. and Munk, A. (2010).
#' \emph{The assessment and planning of non-inferiority trials for retention of effect hypotheses - towards a general approach.}
#' arXiv:0912.4169v1
#' @export
#' @keywords power Poisson
power_RET_poisson <- function(rateExp, rateRef, ratePla, Delta, sig.level = NULL, power = NULL, n = NULL, type = 'unrestricted', allocation = c(1/3, 1/3, 1/3)){

  ## TO-DO
  # - add type='restricted'

  if( sum(sapply(list(n, power, sig.level), is.null)) !=  1 )
    stop("Exactly one of 'n', 'power', and 'sig.level' must be NULL.")

  type <- match.arg(type)

  if (missing(rateExp)) stop("'rateExp' is missing.")
  if (missing(rateRef)) stop("'rateRef' is missing.")
  if (missing(ratePla)) stop("'ratePla' is missing.")
  if (missing(Delta)) stop("'Delta' is missing.")

  if (!is.numeric(c(rateExp, rateRef, ratePla)) || any(0 >= rateExp, 0 >= rateRef,  0 >= ratePla)) {
    stop("'rateExp', 'rateRef', and 'ratePla' must not be larger than 0.")
  }

  # Calculate effect size and check if parameters are located in the alternative
  effect <- rateExp - Delta * rateRef + (Delta - 1) * ratePla
  if (effect >= 0) {
    stop('Parameter vector is not located in the alternative.')
  }

  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 >= sig.level | sig.level >= 1)) {
    stop("'sig.level' must be numeric in (0, 1).")
  }

  if (!is.null(power) && !is.numeric(power) || any(0 >= power | power >= 1)) {
    stop("'power' must be numeric in (0, 1).")
  }

  if (!is.null(Delta) && !is.numeric(Delta) || (0 >= Delta)) {
    stop("'Delta' must be larger than 0.")
  }

  if (!is.null(n) && !is.numeric(n) || any(6 >= n)) {
    stop("'n' must be larger than 6.")
  }

  allocation_check <- any(allocation >= 1, allocation <= 0, sum(allocation) != 1, length(allocation) != 3)
  if ( !is.numeric(allocation) || allocation_check) {
    stop("'allocation' must not have length 3, sum up to 1, and have only entries between 0 and 1.")
  }

  # initialize note
  note <- NULL

  # Adjust Sample Size Allocation
  w <- allocation
  if (!is.null(n)) {
    nExp <- ceiling(n * w[1])
    nRef <- ceiling(n * w[2])
    nPla <- ceiling(n * w[3])
    if (any(c(nExp, nRef, nPla) / sum(c(nExp, nRef, nPla)) != w)) {
      n <-  sum(c(nExp, nRef, nPla))
      w <- c(nExp, nRef, nPla) / n
      note <- "'allocation' and 'n' have been recalculated."
    }
    else if (sum(c(nExp, nRef, nPla)) != n) {
      n <-  sum(c(nExp, nRef, nPla))
      note <- "'n' have been recalculated."
    }
  }

  # Calculate variances
  VarUnres <- rateExp  / w[1] + Delta^2 * rateRef / w[2] + (1 - Delta)^2 * ratePla / w[3]
  VarRes <- VarUnres

  # Define 'method' for output
  switch(type,
         unrestricted = ( method <- 'Wald-type test with unrestricted variance estimation' ),
         restricted = ( method <- 'Wald-type test with restricted maximum-likelihood variance estimation' )
  )

  # Calculate missing parameter
  if (is.null(n)) {
    switch(type,
           unrestricted =  (n <- ceiling( (qnorm(1-sig.level) + qnorm(power))^2 * VarUnres / effect^2) ),
           restricted = (n <- ceiling( (qnorm(1-sig.level)*sqrt(VarRes)/ sqrt(VarUnres) + qnorm(power))^2 * VarUnres / effect^2) )
    )
    nExp <- round(n * w[1])
    nRef <- round(n * w[2])
    nPla <- round(n * w[3])
    n <- nExp + nRef + nPla
    if (any(c(nExp, nRef, nPla) / n != w)) {
      w <- c(nExp, nRef, nPla) / n
      note <- "'allocation' has been recalculated."
    }
  }
  if (is.null(power)) {
    switch(type,
           unrestricted =  (power <- pnorm(qnorm(sig.level) - sqrt(n) * effect / sqrt(VarUnres), mean = 0, sd = 1) ),
           restricted = (power <- pnorm(qnorm(sig.level) * sqrt(VarRes) / sqrt(VarUnres) - sqrt(n) * effect / sqrt(VarUnres), mean = 0, sd = 1))
    )
  }
  if (is.null(sig.level)) {
    switch(type,
           unrestricted =  (sig.level <- pnorm(qnorm(power) + sqrt(n) * effect / sqrt(VarUnres), mean = 0, sd = 1) ),
           restricted = (sig.level <- pnorm(qnorm(power)*sqrt(VarUnres)/ sqrt(VarRes) + sqrt(n) * effect / sqrt(VarRes), mean = 0, sd = 1))
    )
  }

  structure(list(method = method,
                 rateExp = rateExp,
                 rateRef = rateRef,
                 ratePla = ratePla,
                 n = n,
                 sig.level = sig.level,
                 power = power,
                 type = type,
                 Delta = Delta,
                 allocation = w,
                 nExp = nExp,
                 nRef = nRef,
                 nPla = nPla,
                 note = note),
            class = "power.htest")
}




#' @title Power related calcuations for retention of effect hypothesis with Poisson distributed endpoints
#' @description Compute power, sample size, or significance level for Wald-type test
#' testing the retention of effect hypothesis with Poisson distributed endpoints.
#' @details If the individual group sample sizes, i.e. \code{n*allocation} are not natural number,
#' the parameters \emph{n} and \emph{allocation} will be re-calculated.
#' @param x A numeric vector containing the data for the blinded sample size reestimation
#' @param Delta A numeric value specifying the non-inferiority or superiority margin.
#' @param Delta_alternative A numeric value specifying the effect in the alternative
#' @param Delta_PR A numeric value specifying the margin effect
#' \eqn{(\lambda_{P}-\lambda_{E})/(\lambda_{P}-\lambda_{E})} in the alternative
#' @param sig.level A numeric value specifying the significance level (type I error probability)
#' @param power A numeric value specifying the target power (1 - type II error probability)
#' @param allocation A (non-empty) vector specifying the sample size allocation (nExp/n, nRef/n, nPla/n)
#' @param type A character string determing how the variance for the
#' Wald-type test statistic is estimated
#' @return A list with class "power.htest" containing the following components:
#' \item{n}{The total sample size}
#' \item{power}{A numeric value specifying the target power}
#' \item{Delta}{A numeric value specifying the non-inferiority or superiority margin. }
#' \item{sig.level}{A character string specifying the significance level}
#' \item{type}{A character string indicating what type of Wald-type test will be performed}
#' \item{allocation}{A vector with the sample size allocation (nExp/n, nRef/n, nPla/n)}
#' \item{sig.level}{The significance level (Type I error probability)}
#' \item{nExp}{A numeric value specifying the number of sample in the experimental treatment group}
#' \item{nRef}{A numeric value specifying the number of sample in the reference treatment group}
#' \item{nPla}{A numeric value specifying the number of sample in the placebo treatment group}
#' @examples
#' # Initiate rates and margins
#' rateExp <- 0.29
#' rateRef <- 0.2
#' ratePla <- 0.40
#' Delta <- (ratePla - rateExp) / (ratePla - rateRef)
#' Delta_PR <- ratePla - rateRef
#'
#' rateExp <- rateRef
#' Delta_alternative <- 1
#'
#' sig.level <- 0.025
#' power_level <- 0.8
#'
#' # generate random numbers
#' xExp <- rpois(100, lambda = rateExp)
#' xRef <- rpois(100, lambda = rateRef)
#' xPla <- rpois(100, lambda = ratePla)
#' x <- c(xExp, xRef, xPla)
#' # blindly reestimate sample size
#' reestimatesamplesize_RET_poisson(x = x,
#'                                  Delta = Delta,
#'                                  Delta_alternative = Delta_alternative,
#'                                  Delta_PR = Delta_PR,
#'                                  sig.level = sig.level,
#'                                  power = power_level,
#'                                  type = 'unrestricted',
#'                                  allocation = c(1, 1, 1) / 3)
#' @export
#' @references Mielke, M. and Munk, A. (2010).
#' \emph{The assessment and planning of non-inferiority trials for retention of effect hypotheses - towards a general approach.}
#' arXiv:0912.4169v1
#' @keywords power Poisson
reestimatesamplesize_RET_poisson <- function(x, Delta, Delta_alternative, Delta_PR, sig.level, power, type = 'unrestricted', allocation = c(1/3, 1/3, 1/3)){

  ## TO-DO
  # - add type='restricted'

  type <- match.arg(type)

  if (missing(sig.level)) stop("'sig.level' is missing.")
  if (missing(Delta_PR)) stop("'Delta_PR' is missing.")
  if (missing(Delta_alternative)) stop("'Delta_alternative' is missing.")
  if (missing(Delta)) stop("'Delta' is missing.")

  # Calculate effect size and check if parameters are located in the alternative
  effect <- Delta_PR * ( Delta_alternative - Delta )
#   if (effect >= 0) {
#     stop('Parameter vector is not located in the alternative.')
#   }

  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 >= sig.level | sig.level >= 1)) {
    stop("'sig.level' must be numeric in (0, 1).")
  }

  if (!is.null(power) && !is.numeric(power) || any(0 >= power | power >= 1)) {
    stop("'power' must be numeric in (0, 1).")
  }

  if (!is.null(Delta) && !is.numeric(Delta) || (0 >= Delta)) {
    stop("'Delta' must be larger than 0.")
  }
  if (!is.null(Delta_PR) && !is.numeric(Delta_PR) || (0 >= Delta_PR)) {
    stop("'Delta_PR' must be larger than 0.")
  }
  if (!is.null(Delta_alternative) && !is.numeric(Delta_alternative) || (Delta_alternative <= 0)) {
    stop("'Delta_alternative' must be larger than 0.")
  }

  allocation_check <- any(allocation >= 1, allocation <= 0, sum(allocation) != 1, length(allocation) != 3)
  if ( !is.numeric(allocation) || allocation_check) {
    stop("'allocation' must not have length 3, sum up to 1, and have only entries between 0 and 1.")
  }

  # overall mean
  overall_mean <- mean(x)

  wExp <- allocation[1]
  wRef <- allocation[2]
  wPla <- allocation[3]

  # blinded estimation of rates
  rateExp <- overall_mean - Delta_PR * ((wExp + 1) * (1 - Delta_alternative) + wPla)
  rateRef <- overall_mean - Delta_PR * (wExp * (1 - Delta_alternative) + wPla)
  ratePla <- overall_mean + Delta_PR * (1 - wExp * (1 - Delta_alternative) - wPla)


  # Calculate variances
  VarUnres <- rateExp  / wExp + Delta^2 * rateRef / wRef + (1 - Delta)^2 * ratePla / wPla
  VarRes <- VarUnres

  # Define 'method' for output
  switch(type,
         unrestricted = ( method <- 'Wald-type test with unrestricted variance estimation' ),
         restricted = ( method <- 'Wald-type test with restricted maximum-likelihood variance estimation' )
  )

  # Re-calculate missing parameter
  switch(type,
         unrestricted =  (n <- ceiling( (qnorm(1-sig.level) + qnorm(power))^2 * VarUnres / effect^2) ),
         restricted = (n <- ceiling( (qnorm(1-sig.level)*sqrt(VarRes)/ sqrt(VarUnres) + qnorm(power))^2 * VarUnres / effect^2) )
  )

  if (n < length(x)) {
    n <- length(x)
  }

  nExp <- ceiling(n * wExp)
  nRef <- ceiling(n * wRef)
  nPla <- ceiling(n * wPla)
  n <- nExp + nRef + nPla

  #  Check if recalculated sample size fits the input allocation
  note <- NULL
  if (any(c(nExp, nRef, nPla) / n != allocation)) {
    allocation <- c(nExp, nRef, nPla) / n
    note <- "'allocation' has been recalculated."
  }

  structure(list(method = method,
                 n = n,
                 sig.level = sig.level,
                 power = power,
                 type = type,
                 Delta = Delta,
                 Delta_alternative = Delta_alternative,
                 Delta_PR = Delta_PR,
                 allocation = allocation,
                 nExp = nExp,
                 nRef = nRef,
                 nPla = nPla,
                 note = note),
            class = "power.htest")
}


