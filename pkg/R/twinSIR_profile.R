################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### profile-method for class "twinSIR" to calculate the profile log-likelihood
### (normalized) as well as profile likelihood based confidence intervals
###
### Copyright (C) 2009 Michael Hoehle
### $Revision$
### $Date$
################################################################################


######################################################################
# Function to compute likelihood based confidence interval, basically
# the two solutions to
#   f(\theta) = l(\theta)-l(\hat{theta)) + 1/2 dchisq(1-alpha,df=1)=0
# are found.
#
#
# Parameters:
#  logliktilde - normalized likelihood function(theta, ...)
#  theta.hat - the MLE
#  lower - search interval [lower,theta.hat] for f=0
#  upper - search interval [theta.hat,upper] for f=0
#  alpha - confidence level (see Equation 2.6 in Pawitan (2003)
#  ... - additional arguments passed to function logliktilde
######################################################################

likelihood.ci <- function (logliktilde, theta.hat, lower, upper,
    alpha = 0.05, ...)
{
  # Highest Likelihood intervall -- target function
  f <- function(theta, ...) { 
    logliktilde(theta, ...) + 1/2*qchisq(1-alpha, df=1)
  }
  # Compute upper and lower boundary numerically
  hl.lower <- uniroot(f, interval = c(lower, theta.hat), ...)$root
  hl.upper <- uniroot(f, interval = c(theta.hat, upper), ...)$root
  return(c(hl.lower,hl.upper))
}


######################################################################
# Function to compute estimated and profile likelihood based
# confidence intervals. Heavy computations might be necessary!
#
#Params:
# fitted - output from a fit with twinSIR
# profile - list with 4D vector as entries - format:
#               c(index, lower, upper, grid size)
#           where index is the index in the coef vector
#                 lower and upper are the parameter limits (can be NA)
#                 grid size is the grid size of the equally spaced grid
#                 between lower and upper (can be 0)
# alpha - (1-alpha)% profile likelihood CIs are computed.
#         If alpha <= 0 then no CIs are computed
# control - control object to use for optim in the profile loglik computations
#
# Returns:
#  list with profile loglikelihood evaluations on the grid
#  and highest likelihood and wald confidence intervals
######################################################################

profile.twinSIR <- function (fitted, profile, alpha = 0.05,
    control = list(fnscale = -1, factr = 1e1, maxit = 100), ...)
{
  ## Check that input is ok
  profile <- as.list(profile)
  if (length(profile) == 0L) {
    stop("nothing to do")
  }
  lapply(profile, function(one) {
    if (length(one) != 4L) {
      stop("each profile entry has to be of form ",
           "'c(index, lower, upper, grid size)'")
    }})
  if (is.null(fitted[["model"]])) {
    stop("'fitted' must contain the model component")
  }
  
  px <- ncol(fitted$model$X)
  pz <- ncol(fitted$model$Z)
  
  ## Control of the optim procedure
  if (is.null(control[["fnscale",exact=TRUE]])) { control$fnscale <- -1 }
  if (is.null(control[["factr",exact=TRUE]])) { control$factr <- 1e1 }
  if (is.null(control[["maxit",exact=TRUE]])) { control$maxit <- 100 }

  
  ## Estimated normalized likelihood function
  ltildeestim <- function(thetai,i) {
    theta <- theta.ml
    theta[i] <- thetai
    with(fitted$model,
      .loglik(theta, X=X, Z=Z, survs=survs, weights=weights)) - loglik.theta.ml
  }

  ## Profile normalized likelihood function
  ltildeprofile <- function(thetai,i)
  {
    emptyTheta <- rep(0, length(theta.ml))
      
    # Likelihood l(theta_{-i}) = l(theta_i, theta_i)
    ltildethetaminusi <- function(thetaminusi) {
      theta <- emptyTheta
      theta[-i] <- thetaminusi
      theta[i] <- thetai
      with(fitted$model,
       .loglik(theta, X=X, Z=Z, survs=survs, weights=weights)) - loglik.theta.ml
    }
    # Score function of all params except thetaminusi
    stildethetaminusi <- function(thetaminusi) {
      theta <- emptyTheta
      theta[-i] <- thetaminusi
      theta[i] <- thetai
      with(fitted$model,
        .score(theta, X=X, Z=Z, survs=survs, weights=weights))[-i]
    }
      
    # Call optim using L-BFGS-B. For harder constrains we need constr.Optim
    lower <- if (fitted$method == "L-BFGS-B") {
               c(rep(0,px),rep(-Inf,pz))[-i]
             } else {
               -Inf
             }
    upper <- if (fitted$method == "L-BFGS-B") {
               c(rep(Inf,px),rep(Inf,pz))[-i]
             } else {
               Inf
             }
    resOthers <- tryCatch(with(fitted$model,
            optim(theta.ml[-i], fn = ltildethetaminusi, gr = stildethetaminusi,
                  method = fitted$method, control = control,
                  lower = lower, upper = upper)),
            warning = function(w) print(w), error = function(e) list(value=NA))
    resOthers$value
  }
  
  ## Initialize
  theta.ml <- coef(fitted)
  loglik.theta.ml <- c(logLik(fitted))
  se <- sqrt(diag(vcov(fitted)))
  resProfile <- list()
  
  ## Perform profile computations for all requested parameters
  cat("Evaluating the profile logliks on a grid...\n")
  for (i in 1:length(profile))
  {
    cat("i= ",i,"/",length(profile),"\n")
    #Index of the parameter in the theta vector
    idx <- profile[[i]][1]
    #If no borders are given use those from wald intervals (unconstrained)
    if (is.na(profile[[i]][2])) profile[[i]][2] <- theta.ml[idx] - 3*se[idx]
    if (is.na(profile[[i]][3])) profile[[i]][3] <- theta.ml[idx] + 3*se[idx]
    #Evaluate profile loglik on a grid (if requested)
    if (profile[[i]][4] > 0) {
      thetai.grid <- seq(profile[[i]][2],profile[[i]][3],length=profile[[i]][4])
      resProfile[[i]] <- matrix(NA, nrow = length(thetai.grid), ncol = 4L,
        dimnames = list(NULL, c("grid","profile","estimated","wald")))
      
      for (j in 1:length(thetai.grid)) {
        cat("\tj= ",j,"/",length(thetai.grid),"\n")
        resProfile[[i]][j,] <- c(thetai.grid[j],
           ltildeprofile(thetai.grid[j],idx),
           ltildeestim(thetai.grid[j],idx),
#9 June 2009: Bug discovered by L. Held. as part of paper revision. C.f. Pawitan p.63
           - 1/2*(1/se[idx]^2)*(thetai.grid[j] - theta.ml[idx])^2)
      }
    }
  }
#9 June 2009. This did not work.
#  names(resProfile) <- names(theta.ml)[sapply(profile, function(x) x[4L]) > 0]
   names(resProfile) <- names(theta.ml)[sapply(profile, function(x) x[1L])]
  
  ## Profile likelihood intervals
  lower <- if (fitted$method == "L-BFGS-B") {
             c(rep(0,px),rep(-Inf,pz))
           } else {
             -Inf
           }
  ciProfile <- matrix(NA, nrow = length(profile), ncol = 5L,
    dimnames = list(NULL, c("idx","hl.low","hl.up","wald.low","wald.up")))
  if (alpha > 0) {
    cat("Computing likelihood ratio intervals...\n")
    for (i in 1:length(profile))
    {
      cat(i,"/", length(profile),"\n") 
      #Index of the parameter in the theta vector
      idx <- profile[[i]][1]
      #Compute highest likelihood intervals
      ci.hl <- tryCatch(
        likelihood.ci(ltildeprofile, theta.hat = theta.ml[idx],
                      lower = max(lower[idx], theta.ml[idx]-5*se[idx]),
                      upper = theta.ml[idx]+5*se[idx], alpha = alpha, i = idx),
        warning = function(w) print(w),
        error = function(e) rep(NA,2))
      #Wald intervals based on expected fisher information
      ci.wald <- theta.ml[idx] + c(-1,1) * qnorm(1-alpha/2) * se[idx]
      ciProfile[i,] <- c(idx, ci.hl, ci.wald)
    }
    rownames(ciProfile) <- names(theta.ml)[ciProfile[,1]]
  }
  
  return(list(lp=resProfile, ci.hl=ciProfile, profileObj=profile))
}