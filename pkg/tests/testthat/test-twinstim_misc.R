## miscellaneous regression tests for twinstim()
data("imdepi")

## let some districts have no population
imdepi0 <- imdepi
imdepi0$stgrid$popdensity[startsWith(as.character(imdepi$stgrid$tile), "01")] <- 0

## automatic start value is robust against -Inf offset
fit0 <- twinstim(endemic = ~offset(log(popdensity)) + I(start/365),
                 data = imdepi0, model = TRUE)
## beta0 was initialized at Inf in surveillance <= 1.22.1
stopifnot(fit0$converged, is.finite(logLik(fit0)))
