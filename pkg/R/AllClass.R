# -------------  class sts  ----------------------------------------

.sts <- setClass(
    "sts",
    slots = c(
        epoch = "numeric",  # this slot was called "week" in surveillance < 1.3
        freq = "numeric",
        start = "numeric",
        observed = "matrix",
        state = "matrix",
        alarm = "matrix",
        upperbound = "matrix",
        neighbourhood = "matrix",
        populationFrac = "matrix",
        map = "SpatialPolygons",
        control = "list",
        ## New slots added in version 1.1-2 to handle proportion time series:
        epochAsDate = "logical",
        multinomialTS = "logical"
    ),
    prototype = list(
        start = c(2000, 1), freq = 52,  # historical defaults
        epochAsDate = FALSE, multinomialTS = FALSE
    ),
    validity = function (object) {
        dimObserved <- dim(object@observed)
        namesObserved <- colnames(object@observed)
        errors <- c(
            if (!isScalar(object@freq) || object@freq <= 0)
                "'freq' must be a single positive number",
            if (length(object@start) != 2)
                "'start' must be of length two: (year, week/month/idx)",
            if (!is.numeric(object@observed))
                "'observed' must be a numeric matrix",
            ## check consistency of slot dimensions wrt dim(observed):
            if (length(object@epoch) != dimObserved[1L])
                "'epoch' must be of length 'nrow(observed)'",
            if (!identical(dim(object@state), dimObserved))
                "'state' must have the same dimensions as 'observed'",
            if (!identical(dim(object@alarm), dimObserved))
                "'alarm' must have the same dimensions as 'observed'",
            if (!identical(dim(object@upperbound), dimObserved))
                "'upperbound' must have the same dimensions as 'observed'",
            if (!identical(dim(object@neighbourhood), dimObserved[c(2L,2L)]))
                "'neighbourhood' must be a square matrix of size 'ncol(observed)'",
            if (!identical(dim(object@populationFrac), dimObserved))
                "'populationFrac' must have the same dimensions as 'observed'",
            ## disallow NULL colnames in *multivariate* "sts" objects
            if (dimObserved[2L] > 1 && is.null(namesObserved))
                "units must be named (set 'colnames(observed)')",
            ## FIXME: should we generally disallow NULL colnames?
            ## NOTE: aggregate(by="unit") previously (<= 1.15.0) had no colnames
            ## if a map is provided, it must cover all colnames(observed):
            if (length(object@map) > 0 && # i.e., not the empty prototype
                !all(namesObserved %in% row.names(object@map)))
                paste("'map' is incomplete.",
                      "  Ensure that 'all(colnames(observed) %in% row.names(map))'.",
                      sep = "\n"),
            ## check booleans
            if (length(object@epochAsDate) != 1 || is.na(object@epochAsDate))
                "'epochAsDate' must be either TRUE or FALSE",
            ## FIXME: we should enforce epoch[1L] to correspond to start
            ## if (!object@epochAsDate && object@epoch[1L] != 1)
            ##     "'epoch' must be an integer sequence starting at 1",
            if (length(object@multinomialTS) != 1 || is.na(object@multinomialTS))
                "'multinomialTS' must be either TRUE or FALSE"
        )
        ## detect mismatch in column names between different slots
        if (dimObserved[2L] > 1 && !is.null(namesObserved)) {
            slots_dn <- c("state", "alarm", "upperbound", "populationFrac", "neighbourhood")
            errors_dn <- lapply(slots_dn, function (name) {
                cn <- colnames(slot(object, name))
                if (!is.null(cn) && !identical(cn, namesObserved))
                    paste0("'colnames(", name, ")' differ from 'colnames(observed)'")
            })
            errors <- c(errors, unlist(errors_dn))
        }
        if (length(errors) > 0) errors else TRUE
    }
)


######################################################################
# Definition of the stsBP class for backprojections.
######################################################################

setClass("stsBP",
         slots = list(
             ci = "array",
             lambda = "array"
         ),
         contains = "sts")


######################################################################
# Definition of the stsNC class for nowcasts.
######################################################################

setClass("stsNC",
         slots = list(
             reportingTriangle = "matrix",
             predPMF = "list",
             pi = "array",
             truth = "sts",
             delayCDF = "list",
             SR = "array"
         ),
         contains = "sts")
