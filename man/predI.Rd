\name{predI}
\alias{predI}
\alias{print.predI}
\title{Prediction interval for future observations }

\description{
Construct a prediction interval (PCI) for future observations from any continuous distribution. Generic method is \code{print}.
}
\usage{
predI(data, dist, s, n, order=NULL, start=NULL, conf=0.95)

\method{print}{predI}(x, \dots)
}

\arguments{
\item{data}{A numeric vector}
\item{dist}{
A character string \code{"name"} naming a distribution for which the corresponding density function \code{dname}, the corresponding distribution function \code{pname} and the corresponding quantile function \code{qname} must be defined.}
\item{s}{
A numeric vector for the order of the next observation. The length of this vector must be equal to 1. }
\item{n}{
A numeric vector for the size of all data.}
\item{order}{
A numeric vector for the moment order(s). The length of this vector must be equal to the number of parameters to estimate. This argument may be omitted(default) for some distributions for which reasonable order are computed.}
\item{start}{
A named list giving the initial values of parameters of the named distribution. This argument may be omitted(default) for some distributions for which reasonable starting values are computed.}
\item{conf}{Confidence level for the test.}
\item{x}{An object of class \code{"predI"}.  }
\item{\dots}{Further argument to be passed to generic function}
}

\details{
The dist argument is assumed to specify the distribution by the probability density function, the commulative distribution function and the quantile function (d, p, q).
By default, best fitting of the data based on maximum likelihood (mle) and moment matching (mme) methods is performed.
once the parameter(s) is(are) estimated, predI computes the prediction interval (PCI) for the future observation.
This function will be called directly in \code{\link{predP}}.
}

\value{
\code{predI} returns an object of class \code{"predI"}, a list with the following components:
\item{ interval }{ the prediction interval. }
\item{ lower }{ the lower bound of the interval. }
\item{ upper }{ the upper bound of the interval. }
\item{ distname }{ the name of the distribution. }
\item{ r }{ the length of the data. }
\item{ s }{ the order of the next observation. }
\item{ n }{ the length of all the data. }
\item{ parameters }{ the parameter estimate. }

Generic function:
\describe{
\item{\code{print}}{
    The print of a \code{"predI"} object shows few traces about the parameters and the prediction interval.
}}
}
\seealso{ \code{\link{bestfit}}, \code{\link{predP}}.}
\references{
Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.

H. M. Barakat, Magdy E. El-Adll, Amany E. Aly (2014), \emph{Prediction intervals of future observations for a sample random size from any continuous distribution}.
Mathematics and Computers in Simulation, volume 97, 1-13.
}
\author{
H. M. Barakat, O. M. Khaled and Hadeer A. Ghonem.
}

\examples{
# (1) prediction interval for the next observations based on normal distribution
#
set.seed(123)
x1 <- rnorm(15, 2, 4)
predI(x1, "norm", 16, 25)


# (2) prediction interval for the next observations based on weibull distribution
#
library(actuar)
set.seed(123)
x2 <- rweibull(16 , 2 , 3)
predI(x2, "weibull", 20, 20 )

}
