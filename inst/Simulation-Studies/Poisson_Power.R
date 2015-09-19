# library(ThreeArmedTrial)
rm(list = ls())
rateExp <- 0.29
rateRef <- 0.2
ratePla <- 0.40
Delta <- (ratePla - rateExp) / (ratePla - rateRef)
Delta_PR <- ratePla - rateRef

rateExp <- rateRef
Delta_alternative <- (ratePla - rateExp) / (ratePla - rateRef)

sig.level <- 0.025
power_level <- 0.8

number_simlations <- 100000


# Allocation 1:1:1
power_RET_poisson(rateExp = rateExp,
                  rateRef = rateRef,
                  ratePla = ratePla,
                  Delta = Delta,
                  sig.level = sig.level,
                  power = power_level,
                  type = 'unrestricted',
                  allocation = c(1, 1, 1) / 3)

simulate_test_RET_poisson(rateExp = rateExp,
                          rateRef = rateRef,
                          ratePla = ratePla,
                          nExp = 331,
                          nRef = 331,
                          nPla = 331,
                          nSimulations = 500000,
                          Delta = Delta,
                          sig.level = sig.level)


# Allocation 2:1:1
power_RET_poisson(rateExp = rateExp,
                  rateRef = rateRef,
                  ratePla = ratePla,
                  Delta = Delta,
                  sig.level = sig.level,
                  power = power_level,
                  type = 'unrestricted',
                  allocation = c(2, 1, 1) / 4)

simulate_test_RET_poisson(rateExp = rateExp,
                          rateRef = rateRef,
                          ratePla = ratePla,
                          nExp = 468,
                          nRef = 234,
                          nPla = 234,
                          nSimulations = 500000,
                          Delta = Delta,
                          sig.level = sig.level)

# Allocation 2:2:1
power_RET_poisson(rateExp = rateExp,
                  rateRef = rateRef,
                  ratePla = ratePla,
                  Delta = Delta,
                  sig.level = sig.level,
                  power = power_level,
                  type = 'unrestricted',
                  allocation = c(2, 2, 1) / 5)

simulate_test_RET_poisson(rateExp = rateExp,
                          rateRef = rateRef,
                          ratePla = ratePla,
                          nExp = 410,
                          nRef = 410,
                          nPla = 205,
                          nSimulations = 500000,
                          Delta = Delta,
                          sig.level = sig.level)

# Allocation 3:2:1
power_RET_poisson(rateExp = rateExp,
                  rateRef = rateRef,
                  ratePla = ratePla,
                  Delta = Delta,
                  sig.level = sig.level,
                  power = power_level,
                  type = 'unrestricted',
                  allocation = c(3, 2, 1) / 6)

simulate_test_RET_poisson(rateExp = rateExp,
                          rateRef = rateRef,
                          ratePla = ratePla,
                          nExp = 518,
                          nRef = 345,
                          nPla = 172,
                          nSimulations = 500000,
                          Delta = Delta,
                          sig.level = sig.level)
