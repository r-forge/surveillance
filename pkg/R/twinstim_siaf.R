################################################################################
### Spatial interaction functions for twinstim's epidemic component.
### Specific implementations are in separate files (e.g.: Gaussian, power law).
###
### Copyright (C) 2009-2015,2017 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at https://www.R-project.org/Licenses/.
################################################################################



#####################
### "Constructor" ###
#####################

siaf <- function (f, F, Fcircle, effRange, deriv, Deriv, simulate, npars,
                  validpars = NULL)
{
    npars <- as.integer(npars)
    if (length(npars) != 1 || npars < 0L) {
        stop("'siaf$npars' must be a single nonnegative number")
    }
    f <- .checknargs3(f, "siaf$f")
    F <- if (missing(F) || is.null(F)) siaf.fallback.F else {
        F <- match.fun(F)
        if (length(formals(F)) < 4L)
            stop("siaf$F() must accept >=4 arguments ",
                 "(polydomain, f, pars, type)")
        F
    }
    haspars <- npars > 0L
    if (!haspars || missing(deriv)) deriv <- NULL
    if (!is.null(deriv)) deriv <- .checknargs3(deriv, "siaf$deriv")
    if (missing(effRange)) effRange <- NULL
    if (missing(Fcircle) || is.null(Fcircle)) {
        Fcircle <- NULL
        if (!is.null(effRange)) {
            message("'siaf$effRange' only works in conjunction with 'siaf$Fcircle'")
            effRange <- NULL
        }
    }
    if (!is.null(Fcircle)) Fcircle <- .checknargs3(Fcircle, "siaf$Fcircle")
    if (!is.null(effRange)) {
        effRange <- match.fun(effRange)
        if (length(formals(effRange)) < 1L) {
            stop("the 'siaf$effRange' function must accept a parameter vector")
        }
    }
    Deriv <- if (is.null(deriv)) NULL else if (missing(Deriv) || is.null(Deriv))
        siaf.fallback.Deriv else {
            Deriv <- match.fun(Deriv)
            if (length(formals(Deriv)) < 4L)
                stop("siaf$Deriv() must accept >=4 arguments ",
                     "(polydomain, deriv, pars, type)")
            Deriv
        }
    ## Check if simulation function has proper format
    if (missing(simulate)) simulate <- NULL
    if (!is.null(simulate)) {
        simulate <- .checknargs3(simulate, "siaf$simulate")
        if (length(formals(simulate)) == 3L)
            formals(simulate) <- c(formals(simulate), alist(ub=))
    }
    ## Check if the validpars are of correct form
    validpars <- if (!haspars || is.null(validpars))
        NULL else match.fun(validpars)
    ## Done, return result.
    list(f = f, F = F, Fcircle = Fcircle, effRange = effRange,
         deriv = deriv, Deriv = Deriv,
         simulate = simulate,
         npars = npars, validpars = validpars)
}



##########################################
### Constant spatial interaction/dispersal
##########################################

siaf.constant <- function ()
{
    res <- list(
        ## use explicit quote()ing to prevent notes from codetools::checkUsage
        f = as.function(c(alist(s=, pars=NULL, types=NULL),
                          quote(rep.int(1, length(s)/2))),
        ##<- nrow() would take extra time in standardGeneric()
                        envir = .GlobalEnv),
        ## integration over polydomains (F) is handled specially in twinstim
        Fcircle = as.function(c(alist(r=, pars=NULL, type=NULL),
                                quote(pi*r^2)),
                              envir = .GlobalEnv),
        simulate = as.function(c(alist(n=, pars=NULL, type=NULL, ub=),
                                 quote(runifdisc(n, ub))),
                               envir = getNamespace("surveillance")),
        npars = 0L
    )
    attr(res, "constant") <- TRUE
    res
}



##########################################
### Naive defaults for the siaf primitives
##########################################

## numerical integration of f over a polygonal domain (single "owin" and type)
siaf.fallback.F <- function (polydomain, f, pars, type, method = "SV", ...)
{
    if (identical(method,"SV")) {
        polyCub.SV(polyregion = polydomain, f = f, pars, type,
                   alpha = 0, ...) # since max at origin
    } else {
        polyCub(polyregion = polydomain, f = f, method = method,
                pars, type, ...)
    }
}

## numerical integration of f over a circular domain
getFcircle <- function (siaf, control.F = list()) {
    if (is.null(siaf$Fcircle)) {
        function (r, pars, type) {
            disc <- discpoly(c(0,0), r, npoly = 64, class = "owin")
            do.call(siaf$F, c(alist(disc, siaf$f, pars, type), control.F))
        }
    } else {
        siaf$Fcircle
    }
}

## numerical integration of deriv over a polygonal domain
siaf.fallback.Deriv <- function (polydomain, deriv, pars, type,
                                 method = "SV", ...)
{
    deriv1 <- function (s, paridx)
        deriv(s, pars, type)[,paridx,drop=TRUE]
    intderiv1 <- function (paridx)
        polyCub(polyregion = polydomain, f = deriv1, method = method,
                paridx = paridx, ...)
    vapply(X = seq_along(pars), FUN = intderiv1,
           FUN.VALUE = 0, USE.NAMES = FALSE)
}



####################################
### Simulation via polar coordinates (used, e.g., for siaf.powerlaw)
####################################

## Simulate from an isotropic spatial interaction function
## f_{2D}(s) \propto f(||s||), ||s|| <= ub.
## within a maximum distance 'ub' via polar coordinates and the inverse
## transformation method:
## p_{2D}(r,theta) = r * f_{2D}(x,y) \propto r*f(r)
## => angle theta ~ U(0,2*pi) and sample r according to r*f(r)
siaf.simulatePC <- function (intrfr)    # e.g., intrfr.powerlaw
{
    as.function(c(alist(n=, siafpars=, type=, ub=), substitute({
        ## Note: in simEpidataCS, simulation is always bounded to eps.s and to
        ## the largest extend of W, thus, 'ub' is finite
        stopifnot(is.finite(ub))

        ## Normalizing constant of r*f(r) on [0;ub]
        normconst <- intrfr(ub, siafpars, type)

        ## => cumulative distribution function
        CDF <- function (q) intrfr(q, siafpars, type) / normconst

        ## For inversion sampling, we need the quantile function CDF^-1
        ## However, this is not available in closed form, so we use uniroot
        ## (which requires a finite upper bound)
        QF <- function (p) uniroot(function(q) CDF(q)-p, lower=0, upper=ub)$root

        ## Now sample r as QF(U), where U ~ U(0,1)
        r <- vapply(X=runif(n), FUN=QF, FUN.VALUE=0, USE.NAMES=FALSE)
        ## Check simulation of r via kernel estimate:
        ## plot(density(r, from=0, to=ub)); curve(p(x)/normconst,add=TRUE,col=2)

        ## now rotate each point by a random angle to cover all directions
        theta <- runif(n, 0, 2*pi)
        r * cbind(cos(theta), sin(theta))
    })), envir=parent.frame())
}



################################################
### Check F, Fcircle, deriv, Deriv, and simulate
################################################

checksiaf <- function (siaf, pargrid, type = 1, tolerance = 1e-5,
                       method = "SV", ...)
{
    stopifnot(is.list(siaf), is.numeric(pargrid), !is.na(pargrid),
              length(pargrid) > 0)
    pargrid <- as.matrix(pargrid)
    stopifnot(siaf$npars == ncol(pargrid))

    ## Check 'F'
    if (!is.null(siaf$F)) {
        message("'F' vs. cubature using method = \"", method, "\" ... ", appendLF=FALSE)
        comp.F <- checksiaf.F(siaf$F, siaf$f, pargrid, type=type,
                              method=method, ...)
        message(attr(comp.F, "all.equal") <-
            all.equal(comp.F[,1], comp.F[,2],
                      check.attributes=FALSE, tolerance=tolerance))
    }

    ## Check 'Fcircle'
    if (!is.null(siaf$Fcircle)) {
        message("'Fcircle' vs. cubature using method = \"", method, "\" ... ", appendLF=FALSE)
        comp.Fcircle <- checksiaf.Fcircle(siaf$Fcircle, siaf$f, pargrid,
                                          type=type, method=method, ...)
        message(attr(comp.Fcircle, "all.equal") <-
            all.equal(comp.Fcircle[,1], comp.Fcircle[,2],
                      check.attributes=FALSE, tolerance=tolerance))
    }

    ## Check 'deriv'
    if (!is.null(siaf$deriv)) {
        message("'deriv' vs. numerical derivative ... ", appendLF=FALSE)
        if (requireNamespace("maxLik", quietly=TRUE)) {
            maxRelDiffs.deriv <- checksiaf.deriv(siaf$deriv, siaf$f, pargrid,
                                                 type=type)
            message(attr(maxRelDiffs.deriv, "all.equal") <-
                if (any(maxRelDiffs.deriv > tolerance))
                paste("maxRelDiff =", max(maxRelDiffs.deriv)) else TRUE)
        } else message("Failed: need package", sQuote("maxLik"))
    }

    ## Check 'Deriv'
    if (!is.null(siaf$Deriv)) {
        message("'Deriv' vs. cubature using method = \"", method, "\" ... ", appendLF=FALSE)
        comp.Deriv <- checksiaf.Deriv(siaf$Deriv, siaf$deriv, pargrid,
                                      type=type, method=method, ...)
        if (siaf$npars > 1) message() #= LF
        attr(comp.Deriv, "all.equal") <-
            sapply(seq_len(siaf$npars), function (j) {
                if (siaf$npars > 1) message("\tsiaf parameter ", j, ": ", appendLF=FALSE)
                ae <- all.equal(comp.Deriv[,j], comp.Deriv[,siaf$npars+j],
                                check.attributes=FALSE, tolerance=tolerance)
                message(ae)
                ae
            })
    }

    ## Check 'simulate'
    if (interactive() && !is.null(siaf$simulate)) {
        message("Simulating ... ", appendLF=FALSE)
        checksiaf.simulate(siaf$simulate, siaf$f, pargrid[1,], type=type)
        message("(-> check the plot)")
    }

    ## invisibly return check results
    invisible(mget(c("comp.F", "comp.Fcircle",
                     "maxRelDiffs.deriv", "comp.Deriv"),
                   ifnotfound=list(NULL), inherits=FALSE))
}

checksiaf.F <- function (F, f, pargrid, type=1, method="SV", ...)
{
    res <- t(apply(pargrid, 1, function (pars) {
        given <- F(LETTERR, f, pars, type)
        num <- siaf.fallback.F(polydomain = LETTERR, f = f, pars = pars,
                               type = type, method = method, ...)
        c(given, num)
    }))
    colnames(res) <- c("F", method)
    res
}

checksiaf.Fcircle <- function (Fcircle, f, pargrid, type=1,
                               rs=c(1,5,10,50,100), method="SV", ...)
{
    pargrid <- pargrid[rep(1:nrow(pargrid), each=length(rs)),,drop=FALSE]
    rpargrid <- cbind(rs, pargrid, deparse.level=0)
    res <- t(apply(rpargrid, 1, function (x) {
        disc <- discpoly(c(0,0), x[1L], npoly = 128, class = "owin")
        c(ana = Fcircle(x[1L], x[-1L], type),
          num = siaf.fallback.F(polydomain = disc, f = f, pars = x[-1L],
                                type = type, method = method, ...))
    }))
    res
}

checksiaf.deriv <- function (deriv, f, pargrid, type=1, rmax=100)
{
    rgrid <- seq(-rmax,rmax,length.out=21) / sqrt(2)
    rgrid <- rgrid[rgrid != 0] # some siafs are always 1 at (0,0) (deriv=0)
    sgrid <- cbind(rgrid, rgrid)
    apply(pargrid, 1, function (pars) {
        maxLik::compareDerivatives(f, deriv, t0=pars, s=sgrid,
                                   print=FALSE)$maxRelDiffGrad
        ## Note: numDeriv::grad() would only allow one location s at a time
    })
}

checksiaf.Deriv <- function (Deriv, deriv, pargrid, type=1, method="SV", ...)
{
    res <- t(apply(pargrid, 1, function (pars) {
        given <- Deriv(LETTERR, deriv, pars, type)
        num <- siaf.fallback.Deriv(polydomain = LETTERR, deriv = deriv, pars = pars,
                                   type = type, method = method, ...)
        c(given, num)
    }))
    paridxs <- seq_len(ncol(pargrid))
    colnames(res) <- c(paste("Deriv",paridxs,sep="."),
                       paste(method,paridxs,sep="."))
    res
}

checksiaf.simulate <- function (simulate, f, pars, type=1, B=3000, ub=10,
                                plot=interactive())
{
    ## Simulate B points on the disc with radius 'ub'
    simpoints <- simulate(B, pars, type=type, ub=ub)

    if (plot) {
        ## Graphical check in 2D
        opar <- par(mfrow=c(2,1), mar=c(4,3,2,1)); on.exit(par(opar))
        plot(as.im.function(function(x,y,...) f(cbind(x,y), pars, type),
                            W=discpoly(c(0,0), ub, class="owin")),
             axes=TRUE, main="Simulation from the spatial kernel")
        points(simpoints, cex=0.2)
        kdens <- kde2d(simpoints[,1], simpoints[,2], n=100)
        contour(kdens, add=TRUE, col=2, lwd=2,
                labcex=1.5, vfont=c("sans serif", "bold"))
        ##x11(); image(kdens, add=TRUE)

        ## Graphical check of distance distribution
        truehist(sqrt(rowSums(simpoints^2)), xlab="Distance")
        rfr <- function (r) r*f(cbind(r,0), pars, type)
        rfrnorm <- integrate(rfr, 0, ub)$value
        do.call("curve", list(quote(rfr(x)/rfrnorm), add=TRUE, col=2, lwd=2))
        ##<- use do.call-construct to prevent codetools::checkUsage from noting "x"
    }

    ## invisibly return simulated points
    invisible(simpoints)
}
