% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IntegrateD.R
\name{IntegrateD}
\alias{IntegrateD}
\title{Integrate a D-dimensional function}
\usage{
IntegrateD(f, limits)
}
\arguments{
\item{limits}{The limits of the integration. Should be
a named list starting from the outermost integral to the innermost.
Use `quote` or `substitute` to express the limits that vary}

\item{fun}{The function or a one sided formula to be integrated}
}
\value{
the value of the integral
}
\description{
Used to integrate a function that takes in
d parameters and outputs a single value.
}
\examples{
IntegrateD(function(x)1/x,list(x = c(1, 10))) # log(10)
IntegrateD(~x^2+y^2, list(x = c(0, 1), y = c(0, quote(x))))
IntegrateD(\(x,y,z)1,
        list(x = c(0,1), y = c(0, quote(x)), z = c(0, quote(1+x+y))))

## The normal distribution function
fun <- function(x) 1/(2*pi)^(length(x)/2)*exp(-sum(x^2)/2)
# 1-d
all.equal(fun(0.5), dnorm(0.5))
all.equal(IntegrateD(~fun(x), list(x=c(-Inf,0)))[[1]], pnorm(0))
# 2-d
IntegrateD(~fun(c(x,y)), list(x = c(-Inf, 0), y = c(-Inf, 0))) #0.25
IntegrateD(~fun(c(x,y)), list(x = c(-Inf, 0), y = c(quote(x), 0))) #0.25
#Note that the last example cannot be done using `mvtnorm::pmvnorm`

}
