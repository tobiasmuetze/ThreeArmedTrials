#library(ThreeArmedTrials)

rm(list = ls())
number_simlations <- 200000
pvalue <- numeric(number_simlations)
rateExp <- 0.29
rateRef <- 0.2
ratePla <- 0.40
Delta <- (ratePla - rateExp) / (ratePla - rateRef)
Delta_PR <- ratePla - rateRef

sig.level <- 0.025

sqrt( sig.level*(1-sig.level) / number_simlations )

n <- 1050

# 1:1:1
simulate_test_RET_poisson(rateExp = rateExp,
                          rateRef = rateRef,
                          ratePla = ratePla,
                          nExp = n / 3,
                          nRef = n / 3,
                          nPla = n / 3,
                          nSimulations = number_simlations,
                          Delta = Delta,
                          sig.level = sig.level)

# 2:1:1
simulate_test_RET_poisson(rateExp = rateExp,
                          rateRef = rateRef,
                          ratePla = ratePla,
                          nExp = 2 * n / 4,
                          nRef = round(n / 4),
                          nPla = round(n / 4),
                          nSimulations = number_simlations,
                          Delta = Delta,
                          sig.level = sig.level)

# 2:2:1
simulate_test_RET_poisson(rateExp = rateExp,
                          rateRef = rateRef,
                          ratePla = ratePla,
                          nExp = 2 * n / 5,
                          nRef = 2 * n / 5,
                          nPla = n / 5,
                          nSimulations = number_simlations,
                          Delta = Delta,
                          sig.level = sig.level)

# 3:2:1
simulate_test_RET_poisson(rateExp = rateExp,
                          rateRef = rateRef,
                          ratePla = ratePla,
                          nExp = 3 * n / 6,
                          nRef = 2 * n / 6,
                          nPla = n / 6,
                          nSimulations = number_simlations,
                          Delta = Delta,
                          sig.level = sig.level)

