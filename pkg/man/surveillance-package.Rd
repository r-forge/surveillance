%\RdOpts{stage=build}
% Not yet: in R 2.12.0 - 4.0.5, \RdOpts{} had no effect (PR#18073)

\encoding{latin1}
\name{surveillance-package}
\alias{surveillance-package}
\alias{surveillance}
\docType{package}

\title{\pkg{surveillance}: \packageTitle{surveillance}}

\description{
  The \R package \pkg{surveillance} implements statistical methods for the
  retrospective modeling and prospective monitoring of epidemic phenomena
  in temporal and spatio-temporal contexts.
  Focus is on (routinely collected) public health surveillance data,
  but the methods just as well apply to data from environmetrics,
  econometrics or the social sciences. As many of the monitoring methods
  rely on statistical process control methodology, the package is
  also relevant to quality control and reliability engineering.
} 

\details{
  The package implements many typical outbreak detection procedures such
  as Stroup et al. (1989), Farrington et al. (1996), Rossi et al. (1999),
  Rogerson and Yamada (2001), a Bayesian approach (H\enc{�}{oe}hle, 2007),
  negative binomial CUSUM methods (H\enc{�}{oe}hle and Mazick, 2009), and a
  detector based on generalized likelihood ratios (H\enc{�}{oe}hle
  and Paul, 2008), see \code{\link{wrap.algo}}.
  Also CUSUMs for the prospective change-point detection in binomial,
  beta-binomial and multinomial time series are covered based on
  generalized linear modeling, see \code{\link{categoricalCUSUM}}.
  This includes, e.g., paired comparison Bradley-Terry modeling described
  in H\enc{�}{oe}hle (2010), or paired binary CUSUM
  (\code{\link{pairedbinCUSUM}}) described by Steiner et al. (1999).
  The package contains several real-world datasets, the ability
  to simulate outbreak data, visualize the results of the monitoring in
  temporal, spatial or spatio-temporal fashion. In dealing with time
  series data, the fundamental data structure of the package is the S4
  class \code{\link{sts}} wrapping observations, monitoring results and
  date handling for multivariate time series.
  A recent overview of the available monitoring procedures is
  given by Salmon et al. (2016).

  For the retrospective analysis of epidemic spread, the package
  provides three endemic-epidemic modeling frameworks with
  tools for visualization, likelihood inference, and simulation.
  The function \code{\link{hhh4}} offers inference methods for the
  (multivariate) count time series models of Held et al. (2005), Paul et
  al. (2008), Paul and Held (2011), Held and Paul (2012), and Meyer and
  Held (2014). See \code{vignette("hhh4")} for a general introduction
  and \code{vignette("hhh4_spacetime")} for a discussion and
  illustration of spatial \code{hhh4} models.
  Self-exciting point processes are modeled through endemic-epidemic
  conditional intensity functions.
  \code{\link{twinSIR}} (H\enc{�}{oe}hle, 2009) models the
  susceptible-infectious-recovered (SIR) event history of a 
  fixed population, e.g, epidemics across farms or networks;
  see \code{vignette("twinSIR")} for an illustration.
  \code{\link{twinstim}} (Meyer et al., 2012) fits spatio-temporal point
  process models to point patterns of infective events, e.g.,
  time-stamped geo-referenced surveillance data on infectious disease
  occurrence; see \code{vignette("twinstim")} for an illustration.
  A recent overview of the implemented space-time modeling frameworks
  for epidemic phenomena is given by Meyer et al. (2017).
}

%% Author information is dynamically extracted from the DESCRIPTION file
\author{
  \Sexpr[stage=build]{
    pkgdir <- tools:::Rd_macros_package_dir()  # support R CMD Rd2pdf pkg
    desc <- tools:::.read_description(file.path(pkgdir, "DESCRIPTION"))
    aar <- unname(eval(parse(text=desc["Authors@R"])))
    authors <- aar[grep("aut", aar$role)]
    paste0(format(authors, include = c("given", "family")), collapse = ", ")
  }

  Maintainer: Sebastian Meyer \email{seb.meyer@fau.de}
  %\packageMaintainer{surveillance}  % only R >= 4.3.2 knows about Authors@R
}

%% Dynamically extract contributors from the DESCRIPTION file
%% and persons from inst/THANKS for acknowledgement:
\section{Acknowledgements}{
  Substantial contributions of code by:
  \Sexpr[stage=build]{
    contributors <- aar[grepl("ctb", aar$role) & !sapply(aar$family, is.null)]
    paste0(format(contributors, include = c("given", "family")), collapse = ", ")
  }.

  Furthermore, the authors would like to thank the following people
  for ideas, discussions, testing and feedback:
  \Sexpr[stage=build]{
    find_inst_file <- function (file) { # support R CMD Rd2pdf in source package
      if(dir.exists(file.path(pkgdir, "inst")))
        file.path(pkgdir, "inst", file) else
        file.path(pkgdir, file)
    }
    ## internally, Rd objects are assumed to be in UTF-8, see parse_Rd()
    thanks <- readLines(find_inst_file("THANKS"), encoding = "UTF-8")
    toString(grep("^(#|$)", trimws(thanks), invert=TRUE, value=TRUE))
  }.
}

\references{
  \code{citation(package="surveillance")} gives the two main software
  references for the modeling (Meyer et al., 2017) and the monitoring
  (Salmon et al., 2016) functionalities:

  \Sexpr[stage=build,results=rd]{
    paste0("\\\itemize{\n", paste0("\\\item ", tools::toRd(
      readCitationFile(find_inst_file("CITATION"))  # gives marked UTF-8 strings
    ), collapse = "\n\n"), "\n}")
  }

  Further references are listed in \code{surveillance:::REFERENCES}.
  
  If you use the \pkg{surveillance} package in your own
  work, please do cite the corresponding publications.
}

\seealso{
  \url{https://surveillance.R-forge.R-project.org/}
}

\keyword{ package }

\examples{
## Additional documentation and illustrations of the methods are
## available in the form of package vignettes and demo scripts:
vignette(package = "surveillance")
demo(package = "surveillance")
}
