\name{stsplot}
\docType{methods}
\alias{plot,sts,missing-method}
\alias{stsplot} % for convenience

\title{Plot-Methods for Surveillance Time-Series Objects}
\description{
  This page gives an overview of plot types which can be
  produced from objects of class \code{"sts"}.
}

\usage{
  \S4method{plot}{sts,missing}(x, type = observed ~ time | unit, \dots)
}

\arguments{
  \item{x}{an object of class \code{"\linkS4class{sts}"}.}
  \item{type}{see Details.}
  \item{\dots}{arguments passed to the \code{type}-specific plot
    function.}
}

\details{
  There are various types of plots which can be produced from an
  \code{"sts"} object. The \code{type} argument specifies the desired
  plot as a formula, which defaults to \code{observed ~ time | unit},
  i.e., plot the time series of each unit separately. Arguments to
  specific plot functions can be passed as further arguments (\dots).
  The following list describes the plot variants:
  
  \describe{
    \item{\code{observed ~ time | unit}}{The default type shows
      \code{ncol(x)} plots with each containing the time series of one
      observational unit. The actual plotting per unit is done by the
      function \code{\link{stsplot_time1}}, called sequentially from
      \code{\link{stsplot_time}}.
    }
    
    \item{\code{observed ~ time}}{The observations in \code{x} are
      aggregated over units and the resulting univariate time-series
      is plotted. The plotting is done by the function
      \code{\link{stsplot_time}}, which takes the same arguments as the
      old \code{\link{plot.survRes}} function.}
    
    \item{\code{alarm ~ time}}{Generates a so called alarmplot for a
      multivariate \code{sts} object. For each time point and each
      series it is shown whether there is an alarm. In case of
      hierarchical surveillance the user can pass
      an additional argument \code{lvl}, which is a vector of the
      same length as rows in \code{x} specifying for each time series
      its level.
    }

    \item{\code{observed ~ unit}}{
      produces a map of counts (or incidence) per region aggregated over
      time. See \code{\link{stsplot_space}} for optional arguments,
      details and examples.
    }
    \item{\code{observed ~ 1 | unit}}{old version of the map plot, which
      supports shading regions with an alarm.
      The plotting is done by the function \code{\link{stsplot_spacetime}}.
      Use \code{type=observed~unit} for the new implementation as function
      \code{\link{stsplot_space}} (without alarm support, though).
    }

    \item{\code{observed ~ 1 | unit * time}}{old version for animated
      maps via the \code{\link{stsplot_spacetime}} function.
      Each of the \code{nrow(x)} frames contains the
      number of counts per region for the current row in the
      \code{observed} matrix. It is possible to redirect the output
      into files, e.g. to generate an animated GIF.
      NOTE: the new \code{\link{animate.sts}} method supersedes this
      plot \code{type}!
    }
  }
}

\value{
  \code{NULL} (invisibly).
  The methods are called for their side-effects.
}

\seealso{
  the documentation of the individual plot types
  \code{\link{stsplot_time}}, \code{\link{stsplot_space}},
  \code{\link{stsplot_spacetime}} (obsolete), as well as the
  \code{animate}-method \code{\link{animate.sts}}.
  \code{\link{plot.survRes}} is the old implementation.
}

\examples{
data("ha.sts")
print(ha.sts)

## map of total counts by district
plot(ha.sts, type=observed ~ 1 | unit)
## only show a sub-period total for two selected districts
plot(ha.sts[1:20,1:2], type=observed ~ 1 | unit)


\dontrun{
# space-time animation
plot(aggregate(ha.sts,nfreq=13), type= observed ~ 1 | unit * time)

#Configure a png device printer to save the frames
dev.printer <- list(device=png, extension=".png", width=640, height=480,
                    name=file.path(tempdir(),"berlin"))

#Do the animation (without extra sleeping time between frames)
plot(aggregate(ha.sts,nfreq=13), type = observed ~ 1 | unit * time,
     wait.ms=0, dev.printer=dev.printer)

#Use ImageMagick (you might have to adjust the path to 'convert')
system(paste("convert -delay 50 ",dev.printer$name,
           "*.png ", dev.printer$name, "-animated.gif",sep=""))
}
}

\keyword{ts}
\keyword{spatial}
\keyword{hplot}
\keyword{methods}