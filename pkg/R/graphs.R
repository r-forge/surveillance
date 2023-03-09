################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Functions concerning graphs: neighbourhood order, adjacency matrix
### These are wrappers around functionality from package "spdep" by Roger Bivand
###
### Copyright (C) 2009-2013,2017,2023 Sebastian Meyer
### $Revision$
### $Date$
################################################################################


### Determine the matrix of neighbourhood orders
### given the binary matrix of first-order neighbours.
### Working horse: spdep::nblag()

nbOrder <- function (neighbourhood, maxlag = 1)
{
    if (!requireNamespace("spdep"))
        stop("package ", sQuote("spdep"),
             " is required to determine neighbourhood orders")

    stopifnot(isScalar(maxlag), maxlag > 0)
    checkNeighbourhood(neighbourhood)
    neighbourhood <- neighbourhood == 1           # convert to binary matrix
    nregions <- nrow(neighbourhood)
    maxlag <- as.integer(min(maxlag, nregions-1)) # upper bound of nb order

    if (maxlag == 1L) {
        storage.mode(neighbourhood) <- "integer"
        return(neighbourhood)
    }

    ## convert binary adjacency matrix to spdep's "nb" class, "manually",
    ## mimicing spdep::mat2listw(neighbourhood)$neighbours; there is no mat2nb()
    nb <- structure(
        apply(neighbourhood, 1L, function(x)
            if (length(nbs <- which(x))) nbs else 0L,
            simplify = FALSE),
        class = "nb")

    ## compute higher order neighbours using spdep::nblag()
    nb.lags <- spdep::nblag(nb, maxlag=maxlag)

    ## Side note: fast method to determine neighbours _up to_ specific order:
    ## crossprod(neighbourhood) > 0  # up to second order neighbours (+set diag to 0)
    ## (neighbourhood %*% neighbourhood %*% neighbourhood) > 0  # up to order 3
    ## and so on...

    ## convert to a single matrix
    nbmat <- neighbourhood   # logical first-order matrix
    storage.mode(nbmat) <- "numeric"
    for (lag in 2:maxlag) {
        if (any(spdep::card(nb.lags[[lag]]) > 0L)) { # any neighbours of this order
            nbmat.lag <- spdep::nb2mat(nb.lags[[lag]], style="B",
                                       zero.policy=TRUE)
            nbmat <- nbmat + lag * nbmat.lag
        }
    }
    attr(nbmat, "call") <- NULL
    storage.mode(nbmat) <- "integer"

    ## message about maximum neighbour order by region
    maxlagbyrow <- apply(nbmat, 1, max)
    message("Note: range of maximum neighbour order by region is ",
            paste0(range(maxlagbyrow), collapse="-"),
            if (max(maxlagbyrow) == maxlag) " ('maxlag' reached)")

    ## Done
    nbmat
}


### Derive adjacency structure from a SpatialPolygons object
### Working horse: spdep::poly2nb

poly2adjmat <- function (SpP, ..., zero.policy = TRUE)
{
    if (!requireNamespace("spdep"))
        stop("package ", sQuote("spdep"),
             " is required to derive adjacencies from SpatialPolygons")
    nb <- spdep::poly2nb(SpP, ...)
    adjmat <- spdep::nb2mat(nb, style="B", zero.policy=zero.policy)
    attr(adjmat, "call") <- NULL
    colnames(adjmat) <- rownames(adjmat)
    adjmat
}
