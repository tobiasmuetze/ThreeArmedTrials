rm(list = ls())
rateExp <- 0.29
rateRef <- 0.2
ratePla <- 0.40
Delta <- (ratePla - rateExp) / (ratePla - rateRef)
Delta_PR <- ratePla - rateRef

#rateExp <- rateRef
Delta_alternative <- (ratePla - rateExp) / (ratePla - rateRef)
Delta_alternative <- 1

sig.level <- 0.025
power_level <- 0.8


nExp_IPS <- 166
nRef_IPS <- 166
nPla_IPS <- 166
n_IPS <- nExp_IPS + nRef_IPS + nPla_IPS
wExp <- nExp_IPS / n_IPS
wRef <- nRef_IPS / n_IPS
wPla <- nPla_IPS / n_IPS

number_simlations <- 250000
p_value <- n_RSS <- VarUnres <- numeric(number_simlations)


for (i in 1:number_simlations) {

  xExp_IPS <- rpois(nExp_IPS, lambda = rateExp)
  xRef_IPS <- rpois(nRef_IPS, lambda = rateRef)
  xPla_IPS <- rpois(nPla_IPS, lambda = ratePla)

  x <- c(xExp_IPS, xRef_IPS, xPla_IPS)

  RSS <- reestimatesamplesize_RET_poisson(x = x,
                                          Delta = Delta,
                                          Delta_alternative = Delta_alternative,
                                          Delta_PR = Delta_PR,
                                          sig.level = sig.level,
                                          power = power_level,
                                          type = 'unrestricted',
                                          allocation = c(wExp, wRef, wPla))

  nExp <- RSS$nExp
  nRef <- RSS$nRef
  nPla <- RSS$nPla

  n_RSS[i] <- nExp + nRef + nPla

  xExp_RSS <- rpois(nExp-nExp_IPS, lambda = rateExp)
  xRef_RSS <- rpois(nRef-nRef_IPS, lambda = rateRef)
  xPla_RSS <- rpois(nPla-nPla_IPS, lambda = ratePla)

  xExp <- c(xExp_IPS, xExp_RSS)
  xRef <- c(xRef_IPS, xRef_RSS)
  xPla <- c(xPla_IPS, xPla_RSS)

  test_out <- test_RET_poisson(xExp = xExp,
                               xRef = xRef,
                               xPla = xPla,
                               Delta = Delta)

  p_value[i] <- test_out$p.value
}

mean(p_value <= 0.025)

mean(n_RSS)
min(n_RSS)
hist(n_RSS)

simulate_test_RET_poisson(rateExp = rateExp,
                          rateRef = rateRef,
                          ratePla = ratePla,
                          nExp = 331,
                          nRef = 331,
                          nPla = 331,
                          nSimulations = 500000,
                          Delta = Delta,
                          sig.level = sig.level)
