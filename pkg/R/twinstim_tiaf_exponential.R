################################################################################
### Exponential temporal interaction function g(t) = exp(-alpha*t)
###
### Copyright (C) 2009-2014,2017 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at https://www.R-project.org/Licenses/.
################################################################################


## nTypes: determines the number of parameters of the Exponential kernel.
##    In a multitype epidemic, the different types may share
##    the same temporal interaction function (type-invariant), in which case
##    nTypes=1. Otherwise nTypes should equal the number of event types of the
##    epidemic, in which case every type has its own alpha.

tiaf.exponential <- function (nTypes = 1, validpars = NULL)
{
    nTypes <- as.integer(nTypes)
    stopifnot(length(nTypes) == 1L, nTypes > 0L)

    ## function definitions for nTypes = 1 (length(alpha) == 1)
    g <- function (t, alpha, types) {
        exp(-alpha*t)
    }
    G <- function (t, alpha, types) {
        if (alpha==0) t else -exp(-alpha*t)/alpha
    }
    deriv <- function (t, alpha, types) {
        as.matrix( -t*exp(-alpha*t) )
    }
    Deriv <- function (t, alpha, types) {
        as.matrix( if (alpha==0) -t^2/2 else (t+1/alpha)*exp(-alpha*t)/alpha )
    }

    ## adaptions for nTypes > 1
    if (nTypes > 1) {
        ## time points vector t, length(types) = length(t)
        body(g) <- as.call(append(as.list(body(g)),
                                  quote(alpha <- alpha[types]), after=1))
        body(G) <- quote({
            alpha <- alpha[types]
            ifelse (alpha==0, t, -exp(-alpha*t)/alpha)
        })
        body(deriv) <- quote({
            L <- length(t)
            deriv <- matrix(0, L, length(alpha))
            alpha <- alpha[types]
            deriv[cbind(1:L,types)] <- -t*exp(-alpha*t)
            deriv
        })
        body(Deriv) <- quote({
            L <- length(t)
            Deriv <- matrix(0, L, length(alpha))
            alpha <- alpha[types]
            Deriv[cbind(1:L,types)] <-
                ifelse(alpha==0, -t^2/2,
                       (t+1/alpha)*exp(-alpha*t)/alpha)
            Deriv
        })
    }

    ## functions only need the base environment
    environment(g) <- environment(G) <-
        environment(deriv) <- environment(Deriv) <- baseenv()

    ## return the kernel specification
    list(g=g, G=G, deriv=deriv, Deriv=Deriv, npars=nTypes, validpars=validpars)
}
