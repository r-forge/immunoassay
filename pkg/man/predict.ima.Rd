\name{predict.ima}
\alias{predict.ima}
\title{Predict method for "ima" class objects.}
\description{Use this function to obtain predictions from fit data stored in the \code{ima}
    class objects.}
\usage{
    \method{predict}{ima}(object, analyte=1, newdata=NULL, ...)
}
\arguments{
  \item{object}{An object of class \code{ima}.}
  \item{analyte}{Integer in range from 1 to the number of analytes. Determines which 
    analyte will be used for making predictions.}
  \item{newdata}{Numeric vector of MFI values to be used for making predictions. If 
    \code{NULL}, the default calibrator values from the data will be used. This is the
    default.}
  \item{\dots}{Currently ignored.}
}
\details{This function makes predictions using inverse sigmoidal function and the equation
    coefficients stored in the class \code{ima} data frame. This is a very simple function
    and the output is limited to the predictions themselves - i.e. no intervals can be
    obtained.}
\value{Numeric vector of predictions.}
\note{This function works only for class \code{ima} objects that have the coefficients 
    stored in them - that is coming from "xPONENT" platform.}
\author{ Michal J. Figurski, PhD \email{mfigrs@gmail.com}}
\seealso{\code{\link{predict.sigfit}}, \code{\link{print.ima}} and \code{\link{summary.ima}}.}
\examples{\dontrun{
run = read.multiplex("your-path-here")
predict(run)
} }
\keyword{ immunoassay }
\keyword{ ima }
\keyword{ predict }
\keyword{ summary }
\keyword{ fitted object }
