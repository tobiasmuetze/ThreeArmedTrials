context('Poisson distribution')

test_that('Hypothesis test - Error handling',
          {
            expect_error(
              test_RET_poisson(xExp = "a", xRef = rpois(10, lambda = 2),
                               xPla = rpois(10, lambda = 3), Delta = 1),
              "Data must be numeric.")
            expect_error(
              test_RET_poisson(xExp = rpois(10, lambda = 2), xRef = rpois(10, lambda = 2),
                               xPla = rpois(10, lambda = 3), Delta = -1),
              "Delta must be postive.")
            expect_error(
              test_RET_poisson(xExp = rpois(10, lambda = 2), xRef = rpois(10, lambda = 2),
                               xPla = rpois(10, lambda = 3), Delta = 'a'),
              "Delta must be postive.")
          }
)


test_that('Power calculation - Error handling',
          {
            expect_error(
              power_RET_poisson(rateExp = 'a',
                                rateRef = 2,
                                ratePla = 3,
                                Delta = 0.7,
                                sig.level = 0.025,
                                power = 0.8,
                                type = 'unrestricted',
                                allocation = c(1, 1, 1) / 3),
              "'rateExp', 'rateRef', and 'ratePla' must not be larger than 0.")
            expect_error(
              power_RET_poisson(rateExp = 2,
                                rateRef = 2,
                                ratePla = 3,
                                Delta = 0.7,
                                sig.level = 0.025,
                                power = 0.8,
                                type = 'unrestricted',
                                allocation = c(1, 1, 1) / 2),
              "'allocation' must not have length 3, sum up to 1, and have only entries between 0 and 1.")
            expect_error(
              power_RET_poisson(rateExp = 2,
                                rateRef = 2,
                                ratePla = 3,
                                Delta = -1,
                                sig.level = 0.025,
                                power = 0.7,
                                type = 'unrestricted',
                                allocation = c(1, 1, 1) / 3),
              "'Delta' must be larger than 0.")
            expect_error(
              power_RET_poisson(rateExp = 2,
                                rateRef = 2,
                                ratePla = 3,
                                Delta = 0.8,
                                sig.level = 0.025,
                                n = 5,
                                type = 'unrestricted',
                                allocation = c(1, 1, 1) / 3),
              "'n' must be larger than 6.")
          }
)

