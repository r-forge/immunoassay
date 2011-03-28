\name{read.kits}
\alias{read.kits}
\title{Function that reads and processes the kit data csv file.}
\description{
    This function loads the ".csv" file with kit information, sets format for columns and
    calculates midpoints for Quality Control sample ranges.
}
\usage{read.kits(path, file)}
\arguments{  
    \item{path}{Character. A path to the file with kits data.}
    \item{file}{Character. Filename of the kits data file.}
}
\details{
    This function has been pre-defined to a very specific format of input data, and this
    format must be adhered to if this function is to be used. For more information on 
    the format, see help for \code{\link{immunoassay.kits}} data frame.
}
\value{
    This function returns kits \code{\link{immunoassay.kits}} object of class \code{data.frame}.
}
\author{Michal J. Figurski, PhD \email{mfigrs@gmail.com}}
\seealso{\code{\link{read.multiplex}}}
\examples{\dontrun{
    immunoassay.kits = read.kits("your-path", "file name")
}}
\keyword{ immunoassay }
\keyword{ read }
\keyword{ csv }
\keyword{ kit data }

