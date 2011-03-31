\name{immunoassay.kits}
\alias{immunoassay.kits}
\docType{data}
\title{Example dataset with kits data for use by function "sigfit".}
\description{This is an exemplary \code{data.frame} with information on kits used in the 
projects in our laboratory. It is a result of calling function \code{\link{read.kits}} on the
example ".csv" file located in the "immunoassay/data" folder of this package.
Such dataset must be provided in order for the function \code{\link{sigfit}} to work. The 
batch-processing function \code{\link{batch}} will load this data automatically if not 
already present in the search path. See more information there.}
\usage{data(immunoassay.kits)}
\format{
  A data frame with 25 observations on the following 17 variables.
  \describe{
    \item{\code{Lot}}{Factor. Kit lot names / numbers.}
    \item{\code{Analyte}}{Factor. Analyte names.}
    \item{\code{Unit}}{Factor. Analyte units.}
    \item{\code{St1}}{Numeric. Concentration for Standard #1.}
    \item{\code{St2}}{Numeric. Concentration for Standard #2.}
    \item{\code{St3}}{Numeric. Concentration for Standard #3.}
    \item{\code{St4}}{Numeric. Concentration for Standard #4.}
    \item{\code{St5}}{Numeric. Concentration for Standard #5.}
    \item{\code{St6}}{Numeric. Concentration for Standard #6.}
    \item{\code{St7}}{Numeric. Concentration for Standard #7.}
    \item{\code{St8}}{Numeric. Concentration for Standard #8.}
    \item{\code{ConA.lo}}{Numeric. Low concentration limit for Control A.}
    \item{\code{ConA.hi}}{Numeric. High concentration limit for Control A.}
    \item{\code{ConB.lo}}{Numeric. Low concentration limit for Control B.}
    \item{\code{ConB.hi}}{Numeric. High concentration limit for Control B.}
    \item{\code{ConA}}{Numeric. Concentration for Control A.}
    \item{\code{ConB}}{Numeric. Concentration for Control B.}
  }
}
\details{
It is important, due to current limitations of the \code{immunoassay} functions, that the names
of standards and controls are the same in the data files as in this kit file, and that this
convention is kept consistently throughout the project.
}
\source{
"InnoBIA Plasma Abeta forms" and "InnoBIA AlzBio3 CSF biomarkers" kit information.
}
\examples{
data(immunoassay.kits)
str(immunoassay.kits)
}
\keyword{ datasets }
\keyword{ immunoassay }
\keyword{ kit data }
