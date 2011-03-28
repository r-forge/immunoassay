\name{plot.sigfit}
\alias{plot.sigfit}
\alias{print.sigfit}
\alias{summary.sigfit}
\title{Basic display methods for objects fitted by the "sigfit" function.}
\description{Use these functions to print out, summarize and make plots of fitted immunoassay 
    run objects of class \code{sigfit}.}
\usage{
    plot.sigfit(x, table = TRUE, type = "fit", norm = "weighted", ...)
    print.sigfit(x)
    summary.sigfit(x)
}
\arguments{
  \item{x}{An object of class \code{sigfit}.}
  \item{table}{Logical, whether to print the coefficients or not, defaults to \code{TRUE}.}
  \item{type}{Character. Either \code{"fit"} - to plot the fit line and the data points; or
        \code{"resid"} to plot the residuals. Defaults to \code{"fit"}}
  \item{norm}{Character. For the plot \code{type="resid"}, this is the normalization type.
        Either \code{"standardized"} for standardized residuals, or \code{"weighted"} for
        weighted residuals. Defaults to \code{"weighted"}.}
  \item{\dots}{Other parameters passed along to the \code{plot()} function.}
}
\details{
    These are the most basic methods used to display sigmoidal fit objects from the \code{sigfit}
    function. Useful utility functions for LaTeX report templates.
}
\value{
    The functions \code{print.sigfit} and \code{plot.sigfit} are used exclusively for their 
    side effects and do not return anything. Function \code{summary.sigfit} returns a list
    of the following elements:
    \item{fit}{An actual summary of the \code{nls} fit.}
    \item{res}{A table of summary results for Standards and QCs.}
    \item{stats}{A table of fit statistics.}
}
\author{Michal J. Figurski, PhD \email{mfigrs@gmail.com}}
\seealso{\code{\link{sigfit}} and \code{\link{predict.sigfit}}. }
\examples{\dontrun{
    run = read.multiplex("your-path-here")
    fit = sigfit(run)
    fit
    summary(fit)
    plot(fit)
}}
\keyword{ immunoassay }
\keyword{ plot }
\keyword{ print }
\keyword{ summary }
\keyword{ fitted object }
