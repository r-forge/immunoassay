\name{check}
\alias{check}
\title{
Function to examine the sigmoidal fit and provide goodness-of-fit measures.
}
\description{
This is a helper function that returns several error measurement quantities that serve to 
help evaluate the logistic 4-PL and 5-PL models. The model to be checked can either come 
directly from the immunoassay data file or be fitted by the \code{\link{sigfit}} function.
}
\usage{
check(x, analyte, ...)
}
\arguments{
  \item{x}{An object of class \code{sigfit} or \code{ima}.}
  \item{analyte}{Character. The name of the analyte to be evaluated in the \code{class "ima"} 
                 data frame.}
  \item{\dots}{Currently ignored.}
}
\details{
This function is used internally by the \code{\link{batch}}, \code{\link{sigfit}}, 
\code{\link{summary.ima}} and \code{\link{summary.sigfit}} functions to evaluate models. 
}
\value{
    \item{St.error }{Mean and median value of error of fit to Standards (calibrators) 
                     in the measurement scale (not MFI scale).}
    \item{QC.error }{Mean and median value of error of fit to Quality Control samples 
                     in the measurement scale.}
    \item{SSE}{Sum of squares of error values for the fit in the measurement scale. 
               This can be heavily affected by selection of calibrators.}
    \item{Sigma}{Residual sum of squares from the fitted model in the measurement scale. This
                 can be heavily affected by weighting.}
    \item{Syx}{The Sy|x error value for the fit - in the measurement scale.}
    \item{r.squared}{R squared of the fit in the measurement scale. For the \code{nls} 
                fits this can be a negative value in certain situations.}
}
\author{Michal J. Figurski, PhD \email{mfigrs@gmail.com}}
\seealso{\code{\link{sigfit}}, \code{\link{predict.sigfit}}}
\examples{
    \dontrun{
    run = read.multiplex("your-path-here")
    check(run)
    
    fit = sigfit(run)
    check(fit)
    }
}
\keyword{ immunoassay }
\keyword{ multiplex }
\keyword{ sigmoidal }
\keyword{ error }

