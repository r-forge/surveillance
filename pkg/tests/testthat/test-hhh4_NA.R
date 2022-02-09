### fitting hhh4() models to time series with missing values

data("influMen")
fluMen <- disProg2sts(influMen)
## set some observations to NA
set.seed(3)
is.na(fluMen@observed) <- sample(length(fluMen@observed), 100)

## compare endemic-only model against NegBin-GLM
form <- addSeason2formula(f = ~ -1 + fe(1, which = c(TRUE, TRUE)), S = c(3, 1))
fitHHH <- hhh4(fluMen,
               list(end = list(f=form), family = "NegBin1", subset = 1:nrow(fluMen)))
fitGLM <- MASS::glm.nb(
    formula = observed ~ -1 + unit + sin(2*pi*t/52):unit + cos(2*pi*t/52):unit +
        I(sin(4*pi*t/52)*unitI) + I(cos(4*pi*t/52)*unitI) +
        I(sin(6*pi*t/52)*unitI) + I(cos(6*pi*t/52)*unitI),
    data = transform(tidy.sts(fluMen), t = epoch - 1, unitI = unit == "influenza"))

expect_equal(logLik(fitHHH), logLik(fitGLM))
expect_equal(fitted(fitHHH)[!terms(fitHHH)$isNA], unname(fitted(fitGLM)))
expect_equivalent(coef(fitHHH)[["overdisp"]], 1/fitGLM$theta)
idxhhh <- c(1:2, 7:10, 3:6)
expect_equivalent(head(fitHHH$coefficients, -1), fitGLM$coefficients[idxhhh])
expect_equivalent(head(fitHHH$se, -1), summary(fitGLM)$coefficients[idxhhh, 2],
                  tolerance = 0.01)
