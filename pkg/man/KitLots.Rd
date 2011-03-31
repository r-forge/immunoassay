\name{KitLots}
\alias{KitLots}
\docType{data}
\title{Example .csv dataset with kits data.}
\description{This is an exemplary .csv file with kit information - to be read into R using
  \code{\link{read.kits}} function.
}
\format{
  The structure of this file is exactly as the \code{immunoassay.kits} data file, except 
  two last columns with calculated average values for "ConA" and "ConB" samples. See 
  documentation for \code{\link{immunoassay.kits}} for more information.
}
\details{
It is important, due to current limitations of the \code{immunoassay} functions, that the names
of standards and controls are the same in the data files as in this .csv file, and that this
convention is kept consistently throughout the project. Use this file to enter new kits data.
}
\source{
"InnoBIA Plasma Abeta forms" and "InnoBIA AlzBio3 CSF biomarkers" kit information.
}
\examples{\dontrun{
read.kits(path="your-path", file="KitLots.csv")
}}
\keyword{ datasets }
\keyword{ immunoassay }
\keyword{ kit data }
\keyword{ csv }
