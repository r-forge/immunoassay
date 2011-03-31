\name{check}
\alias{check}
\title{
Function to examine the sigmoidal fit done by "sigfit".
}
\description{
This is a helper function that returns several error measurement quantities that serve to 
help evaluate models fitted by the \code{\link{sigfit}} function.
}
\usage{
check(x, ...)
}
\arguments{
  \item{x}{An object of class \code{sigfit}.}
  \item{\dots}{Currently ignored.}
}
\details{
This function is used internally by the \code{\link{batch}} and \code{\link{sigfit}} functions
to evaluate models. 
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
    \item{r.squared}{R squared of the fit in the measurement scale. Keep in mind that this
                     is the \code{nls} fit, and therefore, in certain situations, this can be
                     a negative value.}
}
\author{Michal J. Figurski, PhD \email{mfigrs@gmail.com}}
\seealso{\code{\link{sigfit}}, \code{\link{predict.sigfit}}}
\examples{
    \dontrun{
    run = read.luminex("your-path-here")
    fit = sigfit(run)
    check(fit)
    }
}
\keyword{ immunoassay }
\keyword{ multiplex }
\keyword{ sigmoidal }
\keyword{ error }

