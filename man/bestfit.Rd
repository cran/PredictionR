\name{bestfit}
\alias{bestfit}
\title{Best fitting of a distribution to a data}
\description{
Fit of a distribution to a data by two methods: maximum likelihood (mle) and moment matching (mme). Kolmogorov-Smirnov test is used to construct the best fitting.
}
\usage{
bestfit(data, dist , order=NULL, start=NULL, conf=0.95)
}
\arguments{
\item{data}{A numeric vector}
\item{dist}{
A character string \code{"name"} naming a distribution for which the corresponding density function \code{dname}, the corresponding distribution function \code{pname} and the corresponding quantile function \code{qname} must be defined.}
\item{order}{
A numeric vector for the moment order(s). The length of this vector must be equal to the number of parameters to estimate. This argument may be omitted(default) for some distributions for which reasonable order are computed.}
\item{start}{
A named list giving the initial values of parameters of the named distribution. This argument may be omitted(default) for some distributions for which reasonable starting values are computed.}
\item{conf}{Confidence level for the test.}
}

\details{
This function is not intended to be called directly but is internally called in \code{\link{predI}} and \code{\link{predP}}.
It is assumed that the two methods: "mle" and "mme" are applied then Kolmogorov-Smirnov test is used to construct the best fitting. }
\value{
\code{bestfit} returns a list with following components,
\item{ fit }{ the parameter estimates. }
\item{ p.value }{ the pvalue of the Kolmogorov-Smirnov Test. }
}

\seealso{ \code{\link{predI}}, \code{\link{predP}}.}

\references{
Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.
}

\author{
H. M. Barakat, O. M. Khaled and Hadeer A. Ghonem.
}

\examples{

#best fitting of a logistic distribution
#
n=100
x1 <- rlogis(n, 0.5, 0.8)
bestfit(x1, "logis")
bestfit(x1, "logis")$p.value

}
