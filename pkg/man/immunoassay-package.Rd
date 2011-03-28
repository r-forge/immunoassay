\name{immunoassay-package}
\alias{immunoassay-package}
\alias{immunoassay}
\docType{package}
\title{Functions to aid fitting and batch processing of immunoassay data}
\description{The package \code{immunoassay} contains several functions to aid batch processing of data from
immunoassay runs, although currently only the loading of Multiplex data is implemented. 
There are functions to load data: \code{\link{read.multiplex}}, fit a choice of sigmoidal 
models: \code{\link{sigfit}}, and a batch-processing function to process a whole directory of data files, 
validate the results and produce reports: \code{\link{batch}}. There also are associated \code{print}, 
\code{summary} and \code{plot} methods, that can be used in report templates. Function \code{\link{batch}}
can produce text or LaTeX reports, where the latter are made based upon project-specific templates, 
that are separate R/Sweave scripts.
These functions are in early stage of development - they've been tested on data from a single instrument
platform, in a single study environment, and with 2- and 3-analyte runs only. In other projects 
that utilize different plate layout, or with greater number of analytes, or simply with different 
sample-naming conventions, these functions may require major modifications in order to work properly.
Nevertheless, this package proved remarkably useful in our laboratory allowing for fast processing of large
amounts of data (projects of 70 and more runs) and creation of nice reports and summary plots.
}
\details{
\tabular{ll}{
Package: \tab immunoassay\cr
Type: \tab Package\cr
Version: \tab 0.3\cr
Date: \tab 2011-03-28\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
There are three main functions in this package:
\tabular{ll}{
\code{\link{read.multiplex}}: \tab A function to extract data from a multiplex .csv file.\cr
\code{\link{sigfit}}: \tab A function that fits a choice of sigmoidal models to the data.\cr
\code{\link{batch}}: \tab A batch-processing function that uses the above two functions to process folders of data files.\cr
}
The functions in the \code{immunoassay} package define two classes: \code{ima} class for raw immunoassay data,
and \code{sigfit} class for fitted sigmoidal models. For these classess \code{print}, \code{summary} 
and \code{plot} methods exist. See help on individual functions for more details.

Additionally the \code{immunoassay} package relies on several global vectors:
\tabular{ll}{
immunoassay.kits: \tab A data frame containing kit information (concentrations of Standards and QCs
    for every analyte) for the kits used in producing the data to be processed. This is a required
    element - functions will not run without it, however the \code{\link{batch}} function will
    load the kit data automatically if it exists in proper location and has a proper name. For
    more information, see help on \code{\link{read.kits}} and \code{\link{batch}} functions.\cr
immunoassay.coefs: \tab A list created by \code{\link{batch}} function, that contains equation coefficients
    for every fitted model. This list is used by the fitting function to calculate starting values
    for subsequent fits. This element is not required to start (default starting values will be used), 
    but it will be automatically created. For more information, see help on \code{\link{batch}} function.\cr
immunoassay.options: \tab A list of all the model settings for the data files in the project folder - 
    used by the \code{\link{batch}} function. It is not required to start, but will be automatically created
    by this function. If the processing fails, the records in this list can be modified and re-used.
    For more information, see help on \code{\link{batch}} function.\cr
immunoassay.environment: \tab A pointer to the \code{environment} of the \code{\link{batch}} function, for
    use within the \code{Sweave} project report templates. It is automatically created by the 
    \code{\link{batch}} function. It's necesary due to known shortcoming of \code{Sweave} 
    that can work in the global environment only. \cr
    }
}
\author{Michal J. Figurski, PhD \email{mfigrs@gmail.com} of the Biomarker Research Laboratory, University of
Pennsylvania, Philadelphia, PA.}
\seealso{\code{\link{batch}}, \code{\link{read.multiplex}}, \code{\link{read.kits}} and \code{\link{sigfit}}.}
\keyword{ package }
\keyword{ immunoassay }
\keyword{ batch processing }
\keyword{ multiplex }
\keyword{ sigmoidal }
\keyword{ ADNI }
