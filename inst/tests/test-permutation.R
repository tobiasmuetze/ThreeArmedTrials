context('Permutation test')

test_that('Error handling',
          {
            expect_error(permutation.test(xExp = rnorm(10), xRef = rnorm(10), xPla = rnorm(10),Delta = -1,  nPermutations = 100), "Margin 'Delta' must be postive.")
            expect_error(permutation.test(xExp = rnorm(10), xRef = rnorm(10), xPla = rnorm(10),Delta = 1, nPermutations = 0), "Number of permutations 'nPermutations' must be postive whole number.")
            expect_error(permutation.test(xExp = "a", xRef = rnorm(10), xPla = rnorm(10),Delta = 1, nPermutations = 100), "Input data must be numeric.")
          }
)

