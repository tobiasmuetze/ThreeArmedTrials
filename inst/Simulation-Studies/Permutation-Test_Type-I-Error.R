library(ThreeArmedTrials)
library(xtable)

# Poisson
results_poisson <- data.frame()
# 1:1:1
out <- simulation.permutation.test(parameterExp=list(n = 33, lambda = 0.72),distExp='rpois',
                                   parameterRef=list(n = 33, lambda = 0.34), distRef='rpois',
                                   parameterPla=list(n = 33, lambda = 1.5), distPla='rpois',
                                   nSimulations=15000, nPermutations=10000, Delta=78/116, sig.level=0.025)
results_poisson <- rbind(results_poisson, out)

# 3:2:1
out <- simulation.permutation.test(parameterExp=list(n = 51, lambda = 0.72),distExp='rpois',
                                   parameterRef=list(n = 34, lambda = 0.34), distRef='rpois',
                                   parameterPla=list(n = 17, lambda = 1.5), distPla='rpois',
                                   nSimulations=15000, nPermutations=10000, Delta=78/116, sig.level=0.025)
results_poisson <- rbind(results_poisson, out)

# 2:2:1
out <- simulation.permutation.test(parameterExp=list(n = 40, lambda = 0.72),distExp='rpois',
                                   parameterRef=list(n = 40, lambda = 0.34), distRef='rpois',
                                   parameterPla=list(n = 20, lambda = 1.5), distPla='rpois',
                                   nSimulations=15000, nPermutations=10000, Delta=78/116, sig.level=0.025)
results_poisson <- rbind(results_poisson, out)


results_poisson
xtable(results_poisson, digits = c(0, 0, 2, 0, 2, 0, 2, 5))


# Binary
results_binary <- data.frame()
# 1:1:1
out <- simulation.permutation.test(parameterExp = list(n = 33, size = 1, prob = 0.35), distExp='rbinom',
                                   parameterRef = list(n = 33, size = 1, prob = 0.2), distRef='rbinom',
                                   parameterPla = list(n = 33, size = 1, prob = 0.7), distPla='rbinom',
                                   nSimulations = 15000, nPermutations = 10000, Delta = 0.7, sig.level=0.025)
results_binary <- rbind(results_binary, out)

# 3:2:1
out <- simulation.permutation.test(parameterExp = list(n = 51, size = 1, prob = 0.35), distExp='rbinom',
                                   parameterRef = list(n = 34, size = 1, prob = 0.2), distRef='rbinom',
                                   parameterPla = list(n = 17, size = 1, prob = 0.7), distPla='rbinom',
                                   nSimulations = 15000, nPermutations = 10000, Delta = 0.7, sig.level=0.025)
results_binary <- rbind(results_binary, out)

# 2:2:1
out <- simulation.permutation.test(parameterExp = list(n = 40, size = 1, prob = 0.35), distExp='rbinom',
                                   parameterRef = list(n = 40, size = 1, prob = 0.2), distRef='rbinom',
                                   parameterPla = list(n = 20, size = 1, prob = 0.7), distPla='rbinom',
                                   nSimulations = 15000, nPermutations = 10000, Delta = 0.7, sig.level=0.025)
results_binary <- rbind(results_binary, out)


results_binary
xtable(results_binary, digits = c(0, 0, 0, 2, 0, 0, 2, 0, 0, 2, 5))


0.025 + 2*sqrt(0.025*0.975/15000)
