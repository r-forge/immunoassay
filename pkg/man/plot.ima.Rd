\name{plot.ima}
\alias{plot.ima}
\title{Method for plotting "ima" objects.}
\description{Use these functions to make plots of raw (as well as fitted) immunoassay 
  run objects of class \code{ima}.}
\usage{
    \method{plot}{ima}(x, what="mfi", type="precision", analyte=1, ref=0.25, cts.scale=350, ...)
}
\arguments{
  \item{x}{An object of class \code{ima}.}
  \item{what}{Character, of value either \code{"cts"}, \code{"mfi"}, \code{"res"} or 
    \code{"pred"}. Determines what type of data will be plotted: \code{"cts"} - counts, 
    \code{"mfi"} - MFI, \code{"res"} - results loaded from the data file, \code{"pred"}
    - results obtained from fitting by \code{\link{sigfit}} funtion. Defaults to 
    \code{"mfi"}. 
  \item{type}{Character. This can be either \code{"precision"} - the between-replicate 
    precision, or \code{"accuracy"} - displayed only for samples with nominal concentrations.
  \item{analyte}{Integer of range from 1 to the number of analytes. Determines which analyte will be plotted.}
  \item{ref}{Numeric in range from 0 to 1, determines the size of the reference circle for type \code{"cts"} and 
    \code{"accuracy"} plots. The default is 0.25 (25\%). Use \code{ref=0} to suppress plotting of reference circles.}
  \item{cts.scale}{Numeric. This is the scaling factor for COUNTS only. The default is 350, which made the 
    size of the circles small enough to fit within cells - within our project.}
  \item{\dots}{Other graphical parameters passed along to \code{plot} function.}
}
\details{This function is utilized to display data from \code{ima} objects. It has been 
    tested in limited number of project environments and may not work properly in another 
    setting. In particular the naming of samples and their sequence / placement on the plate
    is critical.}
\value{This functions is used for its side effects and returns nothing.}
\author{Michal J. Figurski, PhD \email{mfigrs@gmail.com}}
\seealso{\code{\link{plot.sigfit}}, \code{\link{print.ima}} and \code{\link{summary.ima}}.}
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
