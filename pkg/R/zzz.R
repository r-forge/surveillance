#######################################
### Hook functions for package start-up
#######################################

gpcWarning <- function ()
    .Deprecated(msg = paste(dQuote("gpc.poly"), "methods are deprecated",
                            "in package", sQuote("surveillance")))

## might no longer be needed in the future:
## https://en.wikipedia.org/wiki/General_Polygon_Clipper
gpclibCheck <- function (fatal = TRUE)
{
    gpclibOK <- surveillance.options("gpclib")
    if (!gpclibOK && fatal) {
        message("Note: The gpclib license is accepted by ",
                sQuote("surveillance.options(gpclib=TRUE)"), ".")
        stop("acceptance of the gpclib license is required")
    }
    gpclibOK
}

.onLoad <- function (libname, pkgname)
{
    ## initialize options
    reset.surveillance.options()
}

.onAttach <- function (libname, pkgname)
{
    ## Startup message
    VERSION <- packageVersion(pkgname, lib.loc=libname)
    packageStartupMessage("This is ", pkgname, " ", VERSION, ". ",
                          "For overview type ",
                          sQuote(paste0("help(", pkgname, ")")), ".")

    if (!interactive()) { # particularly for R CMD check
        ## skip long examples and disallow gpclib, unless:
        allExamples <- nzchar(Sys.getenv("_R_SURVEILLANCE_ALL_EXAMPLES_"))
        ## not using surveillance.options() as this would load gpclib already
        .Options$allExamples$value <- .Options$gpclib$value <- allExamples
    }
}



###########################
### Little helper functions
###########################


### determines multiplicities in a matrix (or data frame)
### and returns unique rows with appended column of counts
### using spatstat's multiplicity methods

countunique <- function (x) unique(cbind(x, COUNT = multiplicity(x)))
