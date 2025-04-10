\encoding{latin1}
\name{stsplot_spacetime}
\alias{stsplot_spacetime}
%\alias{stsplot_spacetime-deprecated}
\title{
  Animated Map of Disease Incidence (DEPRECATED)
}

\description{
  \code{stsplot_spacetime} is \emph{deprecated};
  use the \code{\link[=animate.sts]{animate}} method instead.

  For each period (row) or for the overall period of the
  \code{observed} matrix of the \code{"\linkS4class{sts}"} object, a map
  showing the counts by region is produced.
  It is possible to redirect the output into files, e.g., to generate an
  animated GIF.
}

\usage{
stsplot_spacetime(x, type, legend = NULL, opts.col = NULL, labels = TRUE,
                  wait.ms = 250, cex.lab = 0.7, verbose = FALSE,
                  dev.printer = NULL, ...)
}

\arguments{
  \item{x}{
    an object of class \code{"\linkS4class{sts}"}.
  }
  \item{type}{
    a formula (see \code{\link{stsplot}}). For a map aggregated over
    time (no animation), use \code{observed ~ 1 | unit}, otherwise
    \code{observed ~ 1 | unit * time}.
  }
  \item{legend}{
    a \code{list} containing the following items used
    for coloring
    \describe{
      \item{\code{dx}}{position increments in x direction}
      \item{\code{dy}}{position increments in y direction}
      \item{\code{x}}{position in x}
      \item{\code{y}}{position in y}
      \item{\code{once}}{a Boolean; if \code{TRUE} then only shown once}
    }
    If \code{NULL} then a default legend is used.
  }
  \item{opts.col}{
    a \code{list} containing the two elements
    \describe{
      \item{\code{ncolors}}{number of colors to use for plotting}
      \item{\code{use.color}}{a Boolean; if \code{TRUE} then colors
	will be used in the palette, otherwise grayscale}
    }
  }
  \item{labels}{Boolean whether to add labels
  }
  \item{wait.ms}{Number of milliseconds to wait between each plot
  }
  \item{cex.lab}{\code{cex} of the labels
  }
  \item{verbose}{Boolean whether to write out extra information
  }
  \item{dev.printer}{Either \code{NULL} (default), which means that plotting is
    only to the screen, or a list with elements \code{device},
    \code{extension}, \code{width}, \code{height}, and \code{name}
    (with defaults \code{png}, \code{".png"}, \code{640}, \code{480},
    and \code{"Rplot"}, respectively) to \code{\link{dev.print}} the
    plots to files (only works in interactive sessions).
    This option is more or less obsolete since the \pkg{animation}
    package provides better features for output to files.
  }
  \item{\dots}{Extra arguments sent to the plot function.
  }
}

\author{
  Michael H\enc{�}{oe}hle
}

\seealso{
  Other \code{\link{stsplot}} types, and \code{\link{animate.sts}} for
  the new implementation.
}

\examples{
data("ha.sts")
print(ha.sts)

\dontrun{
# map of total counts by district (compare old vs. new implementation)
plot(ha.sts, type = observed ~ 1 | unit) # deprecated
plot(ha.sts, type = observed ~ unit, labels = TRUE)

# space-time animation
plot(aggregate(ha.sts,nfreq=13), type = observed ~ 1 | unit * time)

#print the frames to a png device
#and do the animation without extra sleeping between frames
imgname <- file.path(tempdir(), "berlin")
plot(aggregate(ha.sts,nfreq=13), type = observed ~ 1 | unit * time,
     wait.ms=0, dev.printer=list(name=imgname))

#Use ImageMagick (you might have to adjust the path to 'convert')
system(paste0("convert -delay 50 ", imgname,
              "*.png ", imgname, "-animated.gif"))
}
}

\keyword{hplot}
\keyword{dynamic}
\keyword{spatial}
\keyword{internal}
