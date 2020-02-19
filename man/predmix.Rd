\name{predmix}
\alias{predmix}
\alias{print.predmix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction point and interval for mixture exponential distribution
}
\description{
Construct a prediction point (Predicted point) and a prediction interval (PI) for mixture exponential distribution. Generic method is \code{print}.
}
\usage{
predmix(data, s, n, parameters,  conf=0.95)

\method{print}{predmix}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A numeric vector}
  \item{s}{
A numeric vector for the order of the next observation. The length of this vector must be equal to 1. }
  \item{n}{
A numeric vector for the size of all data.}
  \item{parameters}{
A numeric vector for the parameters of the distribuion of the data.}
\item{conf}{Confidence level for the test.}
\item{x}{An object of class \code{"predmix"}.  }
\item{\dots}{Further argument to be passed to generic function}
}
\details{
The parameters argument is assumed to specify the parameters for mixture of two exponential distributions by the mixing proportion, rate1 and rate2.
}
\value{
\code{predmix} returns an object of class \code{"predmix"}, a list with the following components:
\item{ point }{ the predicted point. }
\item{ interval }{ the prediction interval. }
\item{ lower }{ the lower bound of the interval. }
\item{ upper }{ the upper bound of the interval. }

Generic function:
\describe{
\item{\code{print}}{
    The print of a \code{"predmix"} object shows few traces about the predicted point and interval.
}}
}
\seealso{ \code{\link{predI}}.}
\references{
H. M. Barakat, O. M. Khaled, Hadeer A. Ghonem (2020), \emph{predicting future
lifetime for mixture exponential distribution}.
Comm. Statist.- Sim. and Comp., DOI: 10.1080/03610918.2020.1715434.

}
\author{
H. M. Barakat, O. M. Khaled and Hadeer A. Ghonem.
}

\examples{
#  prediction for the next observation
#
set.seed(123)
library(Renext)
x<-sort(rmixexp2(30, 0.6, 2, 8))
predmix(x[1:25], 26, 30, c(0.6,2,8))
predmix(x[1:25], 29, 30, c(0.6,2,8))

}
