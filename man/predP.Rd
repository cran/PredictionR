\name{predP}
\alias{predP}
\alias{print.predP}
\title{Prediction point for future observations }

\description{
Construct a prediction point for future observations from any continuous distribution. Generic method is \code{print}.
}
\usage{
predP(data, dist, no, n, order=NULL, start=NULL, conf=0.95)

\method{print}{predP}(x, \dots)
}

\arguments{
\item{data}{A numeric vector}
\item{dist}{
A character string \code{"name"} naming a distribution for which the corresponding density function \code{dname}, the corresponding distribution function \code{pname} and the corresponding quantile function \code{qname} must be defined.}
\item{no}{
A numeric vector for the number of the next observations. The length of this vector plus the length of the data must be less than or equal the length of all  data. }
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
By default, best fitting of the data based on maximum likelihood (mle) and moment matching (mme) methods is performed.
once the parameter(s) is(are) estimated, predP computes the prediction point(s) for the future observation(s).
}
\value{
\code{predP} returns an object of class \code{"predP"}, a list with the following components:
\item{ data }{ the new data with new observations. }
\item{ newobs }{ the new observations. }
\item{ ns }{ the rank of the new observations. }
\item{ no }{ the number of the next observations. }
\item{ distname }{ the name of the distribution. }
\item{ ld }{ the length of the data. }
\item{ n }{ the length of all the data. }

Generic function:
\describe{
\item{\code{print}}{
    The print of a \code{"predP"} object shows the prediction point(s) for the future observation(s).
}}
}

\seealso{ \code{\link{bestfit}}, \code{\link{predI}}.}

\references{
Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.

H. M. Barakat, Magdy E. El-Adll, Amany E. Aly (2014), \emph{Prediction intervals of future observations for a sample random size from any continuous distribution}.
Mathematics and Computers in Simulation, volume 97, 1-13.

H. M. Barakat, O. M. Khaled and Hadeer A. ghonem (2018), \emph{Predicting future lifetime based on random number for mixture exponential distribution}.
International conference of mathematics and its applications (ICMA18), April, 2018, Cairo, Egypt.
}
\author{
H. M. Barakat, O. M. Khaled and Hadeer A. Ghonem.
}

\examples{
# prediction point for the next observations based on gamma distribution
#
set.seed(123)
x1 <- rgamma(10, 4, 2)
predP(x1, "gamma", 8, 20)

}

