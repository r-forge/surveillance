\name{sim.pointSource}
\alias{sim.pointSource}
\encoding{latin1}
\title{Simulate Point-Source Epidemics}
\description{
  Simulation of epidemics which were introduced by point sources.
  The basis of this programme is a combination of a Hidden Markov Model
  (to get random timepoints for outbreaks) and a simple model
  (compare \code{\link{sim.seasonalNoise}}) to simulate the baseline.
}
\usage{
sim.pointSource(p = 0.99, r = 0.01, length = 400, A = 1, 
                alpha = 1, beta = 0, phi = 0, frequency = 1, state = NULL, K)
}
\arguments{
  \item{p}{probability to get a new outbreak at time i if there was one at time i-1, default 0.99.}
  \item{r}{probability to get no new outbreak at time i if there was none at time i-1, default 0.01.}
  \item{length}{number of weeks to model, default 400. \code{length} is ignored if \code{state}
    is given. In this case the length of \code{state} is used.}
  \item{A}{amplitude (range of sinus), default = 1.}
  \item{alpha}{parameter to move along the y-axis (negative values not allowed)
    with alpha > = A, default = 1.}
  \item{beta}{regression coefficient, default = 0.}
  \item{phi}{factor to create seasonal moves
    (moves the curve along the x-axis), default = 0.}
  \item{frequency}{factor to determine the oscillation-frequency, default = 1.}
  \item{state}{use a state chain to define the status at this timepoint (outbreak or not).
    If not given a Markov chain is generated by the programme, default NULL.}
  \item{K}{additional weight for an outbreak which influences the distribution
    parameter mu, default = 0.}
}
\value{
  a \code{disProg} (disease progress) object including a list of the
  observed, the state chain and nearly all input parameters.
}
\seealso{\code{\link{sim.seasonalNoise}}}
\author{M. \enc{H�hle}{Hoehle}, A. Riebler, C. Lang}
\examples{
set.seed(123)
disProgObj <- sim.pointSource(p = 0.99, r = 0.5, length = 208,
                              A = 1, alpha = 1, beta = 0, phi = 0,
                              frequency = 1, state = NULL, K = 2)
plot(disProgObj)

## with predefined state chain
state <- rep(c(0,0,0,0,0,0,0,0,1,1), 20)
disProgObj <- sim.pointSource(state = state, K = 1.2)
plot(disProgObj)

## simulate epidemic, send to RKI 1 system, plot, and compute quality values
testSim <- function (..., K = 0, range = 200:400) {
  disProgObj <- sim.pointSource(..., K = K)
  survResults <- algo.call(disProgObj,
    control = list(list(funcName = "rki1", range = range)))
  plot(survResults[[1]], "RKI 1", "Simulation")
  algo.compare(survResults)
}
testSim(K = 2)
testSim(r = 0.5, K = 5)  # larger and more frequent outbreaks
}
\keyword{datagen}
