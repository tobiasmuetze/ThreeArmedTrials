rm(list = ls())
library(ThreeArmedTrials)

parameters <- read.table(file = "/home/tmuetze/Dropbox/R-Packages/ThreeArmedTrials/inst/Simulation-Studies/Permutation-Test_Type-I-Error_Scenarios.txt", header = TRUE, sep = "\t")

parameters_large <- parameters
parameters_large$rate_exp <- parameters_large$rate_exp * 3
parameters_large$rate_ref <- parameters_large$rate_ref * 3
parameters_large$rate_pla <- parameters_large$rate_pla * 3

parameters_final <- rbind(parameters, parameters_large)
parameters_final$rate_ref <- (parameters_final$rate_exp - (1 - parameters_final$Delta) * parameters_final$rate_pla) / parameters_final$Delta


results <- data.frame()

for (i in seq_along(parameters_final[,1])) {

  n_exp <- parameters_final$n_exp[i]
  rate_exp <- parameters_final$rate_exp[i]
  n_ref <- parameters_final$n_ref[i]
  rate_ref <- parameters_final$rate_ref[i]
  n_pla <- parameters_final$n_pla[i]
  rate_pla <- parameters_final$rate_pla[i]
  Delta <- parameters_final$Delta[i]

  out <- simulate_test_RET_permutation(parameterExp = list(n = n_exp, lambda = rate_exp),
                                       distExp = 'rpois',
                                       parameterRef = list(n = n_ref, lambda = rate_ref),
                                       distRef = 'rpois',
                                       parameterPla = list(n = n_pla, lambda = rate_pla),
                                       distPla = 'rpois',
                                       nSimulations = 50000,
                                       nPermutations = 25000,
                                       Delta = Delta,
                                       sig.level = 0.025,
                                       competitor = TRUE)
  results <- rbind(results, out)
}

results

