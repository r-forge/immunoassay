\name{EXAMPLE-project}
\alias{EXAMPLE-project}
\alias{project}
\alias{project.means}
\alias{project.qcplot}
\alias{example.run-report}
\alias{example.project-report}
\title{Exemplary functions and report templates}
\description{
    In addition to general utility and batch processing functions, the \code{immunoassay} 
    library allows for creating nicely formatted reports. The reports are actually a separate,
    customizable \code{Sweave} scripts that can define their own R functions or utilize the 
    vast array of functions already available in R.
    There are two LaTeX report templates included in this library: 
    \code{"example.run-report.R"} and \code{"example.project-report.R"}. Additionally, there are 
    two functions utilized by these templates: \code{project.means} and \code{project.qcplot} function. 
}
\usage{
project.means(x)
project.qcplot(x, samples=c("ConA","ConB"), type="pred", breaks=c(1,length(unique(x$ID))), log=FALSE, ...)
}
\arguments{
  \item{x}{A \code{data.frame} of results from the \code{\link{batch}} function - 
    either \code{immunoassay.environment$l.run} or \code{immunoassay.environment$results}.}
  \item{samples}{Character vector of sample names to be plotted. A new plot will be called 
    for each sample, so make sure to include \code{par(mfrow=)} parameter when using more
    than one sample in this place. This serves the purpose of having multiple plots on a 
    single pdf page. Sample names must match the names in the data. The default is \code{
    c("ConA","ConB")}.}
  \item{type}{Character. This is the prefix of the column name in the data, that is followed
    by the analyte name. Currently this can have the following values: \code{"pred"} - 
    selects predictions from the model (calculated values); \code{"MFI"} - selects MFI 
    values for plotting; or \code{"CTS"} - selects counts (useless).}
  \item{breaks}{Numeric vector of "breaks" for the calculation and display of \code{mean} 
    and \code{sd} bands on the plot. Useful for projects that utilized more than one kit 
    lot number, to separate results from different kits. Must be provided manually - currently
    no automatic selection of kits is implemented.}
  \item{log}{Logical. If \code{TRUE}, the function takes a natural logarithm of y-axis values.
    Defaults to \code{FALSE}.}
}
\details{
    The two functions mentioned above are simple utility functions called by the EXAMPLE 
    report templates. These templates can be found in the "templates" folder of the  
    \code{\link{immunoassay}} package, along with \code{"Sweave.sty"} LaTeX style. 
    These templates are "Sweaved" in the process of batch-processing of immunoassay data
    using the \code{\link{batch}} function. These files should be copied to the main project 
    folder on the user's computer and modified according to user's needs.
}
\author{Michal J. Figurski, PhD \email{mfigrs@gmail.com} of the Biomarker Research Laboratory, University of
Pennsylvania, Philadelphia, PA.}
\seealso{\code{\link{immunoassay}}, \code{\link{read.multiplex}}, \code{\link{read.kits}} and \code{\link{sigfit}}.}
\keyword{ package }
\keyword{ immunoassay }
\keyword{ batch processing }
\keyword{ multiplex }
\keyword{ sigmoidal }


