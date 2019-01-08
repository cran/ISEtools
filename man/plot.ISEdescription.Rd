% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.ISEdescription.R
\name{plot.ISEdescription}
\alias{plot.ISEdescription}
\title{Plot ISE parameter values}
\usage{
\method{plot}{ISEdescription}(x, breaks = 20, ...)
}
\arguments{
\item{x}{ISE description (e.g. object of class ISEdescription)}
\item{breaks}{Approximate number of bins for histograms, defaults to 20}
\item{...}{Other arguments to be passed through to plotting (histogram) functions}
}
\description{
Plots histograms of ISE parameter values a, b, c, sigma, and LOD (alpha, beta or S/N) for the equation y = a + b log(x + c) + error, where error follows a normal distribution with mean 0 and standard deviation sigma.
}
\seealso{
\code{\link{describeISE}}
}
\author{
Peter Dillingham, \email{peter.dillingham@otago.ac.nz}
}