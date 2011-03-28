\name{predict.sigfit}
\alias{predict.sigfit}
\title{Predict method for "sigfit" function.}
\description{
    This function calculates predictions from the \code{sigfit} models in the scale
    of results. Major difference between this function and predict method for \code{nls} models
    is that \code{predict.nls} would make predictions in the scale of MFI.
}
\usage{
    predict.sigfit(x, newdata = NULL, e.fit = FALSE)
}
\arguments{
    \item{x}{An object of class \code{"sigfit"} - fitted sigmoidal model.}
    \item{newdata}{Data for which to make predictions. If omitted (default), the predictions
        will be made for the standards that were used to fit the model. Otherwise, \code{newdata}
        can be an object of class \code{"ima"}, or a regular \code{data.frame} with MFI 
        data in the 1st column and a placeholder for calculated concentrations (predictions)
        in the second column; other columns are ignored. In case the first column in the 
        \code{data.frame} is named in the fashion of "MFI.analyte-name", the 
        "analyte-name" must match the name of the analyte in the model.}
    \item{e.fit}{Logical. If \code{TRUE}, fit prediction errors (in the percent scale) will
        be added to the result. The default is \code{FALSE}.}
}
\details{
    This function uses the inverse-sigmoid formula to make predictions from \code{nls} model
    in the scale of the result, instead of scale of "MFI".
}
\value{
    If provided with an object of class \code{ima} as an argument to \code{newdata} 
    parameter, this function returns an object of the same type, but with the column of
    predictions added. Otherwise, it returns a numeric vector of predictions (for
    \code{e.fit=FALSE}) or a matrix of predictions and \% error values (for \code{
    e.fit=TRUE}).
}
\author{Michal J. Figurski, PhD \email{mfigrs@gmail.com}}
\seealso{\code{\link{sigfit}}}
\examples{\dontrun{
    run = read.multiplex("your-path-here")
    fit = sigfit(run)
    predict(run)
}}
\keyword{ immunoassay }
\keyword{ predict }
\keyword{ fitted object }
