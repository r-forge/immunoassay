\name{plot.ima}
\alias{plot.ima}
\alias{print.ima}
\title{Methods for plotting and printing "ima" objects.}
\description{Use these functions to print out and make plots of raw (as well as fitted) immunoassay 
  run objects of class \code{ima}.}
\usage{
plot.ima(data, type = "cts", analyte = 1, ref = 0.25, cts.scale = 350, ...)
print.ima(data)
}
\arguments{
  \item{data}{An object of class \code{ima}.}
  \item{type}{Character. This can be either \code{"cts"}, \code{"cv"}, or \code{"accuracy"}. 
    The \code{type="cts"} is the default for un-fitted \code{ima} data objects and it plots COUNTS for visual
    inspection. The \code{type="cv"} for the un-fitted \code{ima} object displays the \%CV of the MFIs, 
    and for fitted \code{ima} objects (must contain predictions from the model), displays the \%CV
    of the results. The \code{type="accuracy"} is an extension of capability of this function for fitted 
    immunoassay data objects and it plots accuracy of the predictions for Standards and QCs.}
  \item{analyte}{Integer of range from 1 to the number of analytes. Determines which analyte will be plotted.}
  \item{ref}{Numeric in range from 0 to 1, determines the size of the reference circle for type \code{"cts"} and 
    \code{"accuracy"} plots. The default is 0.25 (25\%). Use \code{ref=0} to suppress plotting of reference circles.}
  \item{cts.scale}{Numeric. This is the scaling factor for COUNTS only. The default is 350, which made the 
    size of the circles small enough to fit within cells - within our project.}
  \item{\dots}{Other graphical parameters passed along to \code{plot} function.}
}
\details{These functions are utilized to display data from \code{ima} objects and additionally provide some functionality
for fitted objects. These functions have been thoroughly tested in only one project environment and may not work
properly in another setting. In particular the naming of samples and their sequence / placement on the plate
is critical.}
\value{These functions are used for their side effects.}
\author{Michal J. Figurski, PhD \email{mfigrs@gmail.com}}
\seealso{\code{\link{plot.sigfit}}, \code{\link{print.sigfit}} and \code{\link{summary.sigfit}}.}
\examples{\dontrun{
run = read.multiplex("your-path-here")
plot(run)
run
} }
\keyword{ immunoassay }
\keyword{ ima }
\keyword{ plot }
\keyword{ print }
\keyword{ fitted object }
