context('Wald-type test for normal distribution: test_RET_normal')

test_that('Error and warning handling',
          {
            expect_error(
              test_RET_normal(xExp = rnorm(10),
                              xRef = rnorm(10),
                              xPla = rnorm(10),
                              Delta = -1),
              "Margin 'Delta' must be postive."
            )
            expect_error(
              test_RET_normal(xExp = rnorm(10),
                              xRef = rnorm(10),
                              xPla = rnorm(10),
                              Delta = "a"),
              "Margin 'Delta' must be a numeric variable."
            )
            expect_error(
              test_RET_normal(xExp = "a",
                              xRef = rnorm(10),
                              xPla = rnorm(10),
                              Delta = 0.8),
              "'xExp', 'xRef', must 'xPla' contain only numeric entries."
            )
            expect_error(
              test_RET_normal(xExp = rnorm(10),
                              xRef = rnorm(10),
                              xPla = rnorm(10),
                              Delta = 0.8,
                              var.equal = "1"),
              "'var.equal' is not a logical variable."
            )
            expect_error(
              test_RET_normal(xExp = rnorm(10),
                              xRef = rnorm(10),
                              xPla = rnorm(10),
                              Delta = 0.8,
                              var.equal = 0),
              "'var.equal' is not a logical variable."
            )
            expect_error(
              test_RET_normal(xExp = rnorm(10),
                              xRef = rnorm(10),
                              xPla = rnorm(10),
                              Delta = 0.8,
                              var.equal = 1),
              "'var.equal' is not a logical variable."
            )
            expect_warning(
              test_RET_normal(xExp = c(NA, rnorm(2)),
                              xRef = rnorm(10),
                              xPla = rnorm(10),
                              Delta = 0.8),
              "'xExp', 'xRef', or 'xPla' contain NA entries which will be removed."
            )
          }
)

