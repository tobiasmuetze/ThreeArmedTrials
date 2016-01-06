rm(list = ls())
library("ThreeArmedTrials")
library("parallel")
library("foreach")
library("doMC")


# ##############################################################
# # Simulations for Poisson Data
# ##############################################################
#
# # Read parameters
# rm(list = ls())
#
# path_parameters <- system.file("Simulation-Studies", "Permutation_Test", "Type-I-Error_Scenarios_Poisson.txt", package = "ThreeArmedTrials")
# parameters <- read.table(file = path_parameters, header = TRUE, sep = "\t")
#
# # Set working directory and output file name
# work_dir <- "./"
# output_filename <- paste0(work_dir, "Results", basename(path_parameters))
#
# # Parallel computation using only one core
# registerDoMC( cores = detectCores() - 1)
#
# results <- foreach(i = seq_along(parameters[,1]), .combine = rbind) %dopar% {
#
#   n_exp <- parameters$n_exp[i]
#   rate_exp <- parameters$rate_exp[i]
#   n_ref <- parameters$n_ref[i]
#   rate_ref <- parameters$rate_ref[i]
#   n_pla <- parameters$n_pla[i]
#   rate_pla <- parameters$rate_pla[i]
#   Delta <- parameters$Delta[i]
#   rate_exp - Delta * rate_ref + (Delta - 1) * rate_pla
#   out <- simulate_test_RET_permutation(parameterExp = list(n = n_exp, lambda = rate_exp),
#                                        distExp = 'rpois',
#                                        parameterRef = list(n = n_ref, lambda = rate_ref),
#                                        distRef = 'rpois',
#                                        parameterPla = list(n = n_pla, lambda = rate_pla),
#                                        distPla = 'rpois',
#                                        nSimulations = 1000,
#                                        nPermutations = 1000,
#                                        Delta = Delta,
#                                        sig.level = 0.025,
#                                        competitor = TRUE)
#   out
#
# }
#
# write.table(results, file = output_filename, sep = "\t", row.names = FALSE, col.names = TRUE)
# ##############################################################
# # END Simulations for Poisson Data
# ##############################################################
#
#
#
# ##############################################################
# # Simulations for negative binomial data - Var = 2 * Exp
# ##############################################################
# rm(list = ls())
#
# # Read parameters
# path_parameters <- system.file("Simulation-Studies", "Permutation_Test",
#                                "Type-I-Error_Scenarios_Poisson.txt",
#                                package = "ThreeArmedTrials")
#
# parameters <- read.table(file = path_parameters, header = TRUE, sep = "\t")
#
# # Set working directory and output file name
# work_dir <- "./"
# output_filename <- paste0(work_dir, "Results-Type-I-Error_Scenarios_NegBin_20.txt")
#
# # Parallel computation using only one core
# registerDoMC( cores = detectCores() - 1)
#
# results <- foreach(i = seq_along(parameters[,1]), .combine = rbind) %dopar% {
#
#   n_exp <- parameters$n_exp[i]
#   rate_exp <- parameters$rate_exp[i]
#   n_ref <- parameters$n_ref[i]
#   rate_ref <- parameters$rate_ref[i]
#   n_pla <- parameters$n_pla[i]
#   rate_pla <- parameters$rate_pla[i]
#   Delta <- parameters$Delta[i]
#   rate_exp - Delta * rate_ref + (Delta - 1) * rate_pla
#   shape <- parameters$shape[i]
#   out <- simulate_test_RET_permutation(parameterExp = list(n = n_exp, mu = rate_exp, size = rate_exp),
#                                        distExp = 'rnbinom',
#                                        parameterRef = list(n = n_ref, mu = rate_ref, size = rate_ref),
#                                        distRef = 'rnbinom',
#                                        parameterPla = list(n = n_pla, mu = rate_pla, size = rate_pla),
#                                        distPla = 'rnbinom',
#                                        nSimulations = 1000,
#                                        nPermutations = 1000,
#                                        Delta = Delta,
#                                        sig.level = 0.025,
#                                        competitor = TRUE)
#   out
#
# }
#
# write.table(results, file = output_filename, sep = "\t", row.names = FALSE, col.names = TRUE)
# ##############################################################
# # END Simulations for negative binomial data - Var = 2 * Exp
# ##############################################################
#
#
#
#
# ##############################################################
# # Simulations for negative binomial data - Var = 3 * Exp
# ##############################################################
# rm(list = ls())
#
# # Read parameters
# path_parameters <- system.file("Simulation-Studies", "Permutation_Test",
#                                "Type-I-Error_Scenarios_Poisson.txt",
#                                package = "ThreeArmedTrials")
#
# parameters <- read.table(file = path_parameters, header = TRUE, sep = "\t")
#
# # Set working directory and output file name
# work_dir <- "./"
# output_filename <- paste0(work_dir, "Results-Type-I-Error_Scenarios_NegBin_30.txt")
#
# # Parallel computation using only one core
# registerDoMC( cores = detectCores() - 1)
#
# results <- foreach(i = seq_along(parameters[,1]), .combine = rbind) %dopar% {
#
#   n_exp <- parameters$n_exp[i]
#   rate_exp <- parameters$rate_exp[i]
#   n_ref <- parameters$n_ref[i]
#   rate_ref <- parameters$rate_ref[i]
#   n_pla <- parameters$n_pla[i]
#   rate_pla <- parameters$rate_pla[i]
#   Delta <- parameters$Delta[i]
#   rate_exp - Delta * rate_ref + (Delta - 1) * rate_pla
#   shape <- parameters$shape[i]
#   out <- simulate_test_RET_permutation(parameterExp = list(n = n_exp, mu = rate_exp, size = rate_exp/2),
#                                        distExp = 'rnbinom',
#                                        parameterRef = list(n = n_ref, mu = rate_ref, size = rate_ref/2),
#                                        distRef = 'rnbinom',
#                                        parameterPla = list(n = n_pla, mu = rate_pla, size = rate_pla/2),
#                                        distPla = 'rnbinom',
#                                        nSimulations = 1000,
#                                        nPermutations = 1000,
#                                        Delta = Delta,
#                                        sig.level = 0.025,
#                                        competitor = TRUE)
#   out
#
# }
#
# write.table(results, file = output_filename, sep = "\t", row.names = FALSE, col.names = TRUE)
# ##############################################################
# # END Simulations for negative binomial data - Var = 3 * Exp
# ##############################################################
#






##############################################################
# Simulations for log-normal data
##############################################################
# rm(list = ls())
#
# # Read parameters
# path_parameters <- system.file("Simulation-Studies",
#                                "Permutation_Test",
#                                "Type-I-Error_Scenarios_Poisson.txt",
#                                package = "ThreeArmedTrials")
# parameters <- read.table(file = path_parameters, header = TRUE, sep = "\t")
#
# # Set working directory and output file name
# work_dir <- "./"
# output_filename <- paste0(work_dir, "Results lognormal - half variance.txt")
#
# # Parallel computation using only one core
# registerDoMC( cores = detectCores() - 1)
#
# results <- foreach(i = seq_along(parameters[,1]), .combine = rbind) %dopar% {
#
#   n_exp <- parameters$n_exp[i]
#   n_ref <- parameters$n_ref[i]
#   n_pla <- parameters$n_pla[i]
#
#   sdlog_exp <- sqrt(log(1 / (2 * parameters$rate_exp[i]) + 1))
#   meanlog_exp <- log(parameters$rate_exp[i]) - sdlog_exp^2/ 2
#
#   sdlog_ref <- sqrt(log(1 / (2 * parameters$rate_ref[i]) + 1))
#   meanlog_ref <- log(parameters$rate_ref[i]) - sdlog_ref^2 / 2
#
#   sdlog_pla <- sqrt(log(1 / (2 * parameters$rate_pla[i]) + 1))
#   meanlog_pla <- log(parameters$rate_pla[i]) - sdlog_pla^2 / 2
#
#   Delta <- parameters$Delta[i]
#
#   out <- simulate_test_RET_permutation(parameterExp = list(n = n_exp, meanlog = meanlog_exp, sdlog = sdlog_exp),
#                                        distExp = 'rlnorm',
#                                        parameterRef = list(n = n_ref, meanlog = meanlog_ref, sdlog = sdlog_ref),
#                                        distRef = 'rlnorm',
#                                        parameterPla = list(n = n_pla, meanlog = meanlog_pla, sdlog = sdlog_pla),
#                                        distPla = 'rlnorm',
#                                        nSimulations = 250,
#                                        nPermutations = 150,
#                                        Delta = Delta,
#                                        sig.level = 0.025,
#                                        competitor = TRUE)
#   out
#
# }
#
# write.table(results, file = output_filename, sep = "\t", row.names = FALSE, col.names = TRUE)
##############################################################
# END Simulations for log-normal data
##############################################################




##############################################################
# Simulations for log-normal data - Var = Exp
##############################################################
# rm(list = ls())
#
# # Read parameters
# path_parameters <- system.file("Simulation-Studies",
#                                "Permutation_Test",
#                                "Type-I-Error_Scenarios_Poisson.txt",
#                                package = "ThreeArmedTrials")
# parameters <- read.table(file = path_parameters, header = TRUE, sep = "\t")
#
# # Set working directory and output file name
# work_dir <- "./"
# output_filename <- paste0(work_dir, "Results-Type-I-Error_Lognorm-VarExp.txt")
#
# # Parallel computation using only one core
# registerDoMC( cores = detectCores() - 1)
#
# results <- foreach(i = seq_along(parameters[,1]), .combine = rbind) %dopar% {
#
#   n_exp <- parameters$n_exp[i]
#   n_ref <- parameters$n_ref[i]
#   n_pla <- parameters$n_pla[i]
#
#   sdlog_exp <- sqrt(log(1 / (parameters$rate_exp[i]) + 1))
#   meanlog_exp <- log(parameters$rate_exp[i]) - sdlog_exp^2/ 2
#
#   sdlog_ref <- sqrt(log(1 / (parameters$rate_ref[i]) + 1))
#   meanlog_ref <- log(parameters$rate_ref[i]) - sdlog_ref^2 / 2
#
#   sdlog_pla <- sqrt(log(1 / (parameters$rate_pla[i]) + 1))
#   meanlog_pla <- log(parameters$rate_pla[i]) - sdlog_pla^2 / 2
#
#   Delta <- parameters$Delta[i]
#
#   out <- simulate_test_RET_permutation(parameterExp = list(n = n_exp, meanlog = meanlog_exp, sdlog = sdlog_exp),
#                                        distExp = 'rlnorm',
#                                        parameterRef = list(n = n_ref, meanlog = meanlog_ref, sdlog = sdlog_ref),
#                                        distRef = 'rlnorm',
#                                        parameterPla = list(n = n_pla, meanlog = meanlog_pla, sdlog = sdlog_pla),
#                                        distPla = 'rlnorm',
#                                        nSimulations = 25000,
#                                        nPermutations = 15000,
#                                        Delta = Delta,
#                                        sig.level = 0.025,
#                                        competitor = TRUE)
#   out
#
# }
#
# write.table(results, file = output_filename, sep = "\t", row.names = FALSE, col.names = TRUE)
##############################################################
# Simulations for log-normal data - Var = Exp
##############################################################




##############################################################
# Simulations for binary data
##############################################################
rm(list = ls())

# Read parameters
path_parameters <- system.file("Simulation-Studies",
                               "Permutation_Test",
                               "Type-I-Error_Scenarios_Bernoulli.txt",
                               package = "ThreeArmedTrials")
parameters <- read.table(file = path_parameters, header = TRUE, sep = "\t")

# Set working directory and output file name
work_dir <- "./"
output_filename <- paste0(work_dir, "Results-Type-I-Error_Bernoulli.txt")

# Parallel computation using only one core
registerDoMC( cores = detectCores() - 1)

results <- foreach(i = seq_along(parameters[,1]), .combine = rbind) %dopar% {

  n_exp <- parameters$n_exp[i]
  n_ref <- parameters$n_ref[i]
  n_pla <- parameters$n_pla[i]

  prob_exp <- parameters$rate_exp[i]
  prob_ref <- parameters$rate_ref[i]
  prob_pla <- parameters$rate_pla[i]

  Delta <- parameters$Delta[i]

  out <- simulate_test_RET_permutation(parameterExp = list(n = n_exp, size = 1, prob = prob_exp),
                                       distExp = 'rbinom',
                                       parameterRef = list(n = n_ref, size = 1, prob = prob_ref),
                                       distRef = 'rbinom',
                                       parameterPla = list(n = n_pla, size = 1, prob = prob_pla),
                                       distPla = 'rbinom',
                                       nSimulations = 14000,
                                       nPermutations = 100,
                                       Delta = Delta,
                                       sig.level = 0.025,
                                       competitor = TRUE)
  out

}

write.table(results, file = output_filename, sep = "\t", row.names = FALSE, col.names = TRUE)
##############################################################
# END Simulations for binary data
##############################################################



# # Non-parallel computation using only one core
# for (i in seq_along(parameters[,1])) {
#
#   n_exp <- parameters$n_exp[i]
#   rate_exp <- parameters$rate_exp[i]
#   n_ref <- parameters$n_ref[i]
#   rate_ref <- parameters$rate_ref[i]
#   n_pla <- parameters$n_pla[i]
#   rate_pla <- parameters$rate_pla[i]
#   Delta <- parameters$Delta[i]
#
#   out <- simulate_test_RET_permutation(parameterExp = list(n = n_exp, lambda = rate_exp),
#                                        distExp = 'rpois',
#                                        parameterRef = list(n = n_ref, lambda = rate_ref),
#                                        distRef = 'rpois',
#                                        parameterPla = list(n = n_pla, lambda = rate_pla),
#                                        distPla = 'rpois',
#                                        nSimulations = 50000,
#                                        nPermutations = 25000,
#                                        Delta = Delta,
#                                        sig.level = 0.025,
#                                        competitor = TRUE)
#   results <- rbind(results, out)
# }
#
