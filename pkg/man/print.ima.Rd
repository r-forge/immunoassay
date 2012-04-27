\name{print.ima}
\alias{print.ima}
\alias{summary.ima}
\title{Methods for printing "ima" objects.}
\description{Use these functions to print out data and summary of immunoassay run objects 
    of class \code{ima}.}
\usage{
    \method{print}{ima}(x, ...)
    \method{summary}{ima}(x, analyte="all", result="res", type="fit")
}
\arguments{
  \item{x}{An object of class \code{ima}.}
  \item{analyte}{Either character or integer in range from 1 to the number of analytes. 
    Determines which analyte will be plotted. If \code{"all"} (the default) - the summary
    for all analytes in sequence is printed. If numeric - only the summary for given analyte
    is displayed.}
  \item{result}{Character. Select which type of result is to be summarized. Use \code{"res"} 
    for results loaded from the run file; or \code{"pred"} for results obtained by fitting 
    the data with the \code{\link{sigfit}} function.}
  \item{type}{Character. This can be either \code{"fit"} or \code{"data"}. For \code{"fit"},
    the fit information is summarized. For \code{"data"} the summary of the results \code{data.frame}
    is displayed.}
  \item{\dots}{Other parameters passed along to \code{print} function.}
}
\details{These functions are utilized to display and summarize data from \code{ima} objects.}
\value{These functions are used for their side effects and return nothing.}
\author{Michal J. Figurski, PhD \email{mfigrs@gmail.com}}
\seealso{\code{\link{plot.sigfit}}, \code{\link{print.sigfit}} and \code{\link{summary.sigfit}}.}
\examples{\dontrun{
run = read.multiplex("your-path-here")
head(run) # Makes use of the plot.ima() function.
summary(run, an=1)
} }
\keyword{ immunoassay }
\keyword{ ima }
\keyword{ print }
\keyword{ summary }
\keyword{ fitted object }
