\name{read.multiplex}
\alias{read.multiplex}
\title{Load raw data from multiplex (x-MAP technology) run file.}
\description{
    This fuction reads a multiplex run .csv file from a given folder and extracts raw data 
    from it (MFI, counts, results and fit information, if available). It returns an object 
    of class \code{ima} which is a specifically formatted \code{data.frame} with immunoassay 
    run parameters stored as attributes.
}
\usage{ read.multiplex(path, file, analytes = "all") }
\arguments{
  \item{path}{Character. The path to the data file. }
  \item{file}{Character. Filename of the file to load. Does not need to contain ".csv" 
    extension - it will be added automatically. }
  \item{analytes}{Either character \code{'analytes="all"'} to load all available analytes, 
    or a numeric vector, to select from the analytes available in the data.}
}
\details{
    This function reads the raw data from multiplex run files: MFIs, COUNTs results, run 
    information. It can also load nominal concentrations of Standards and QCs along with 
    fit information, if available (currently "xPONENT" files only). 
    
    It then creates a specifically formatted data frame with many run parameters stored as
    attributes. It can currently recognize two ".csv" file formats - the "Luminex IS" and
    "xPONENT" output. If your file doesn't load properly, please email it to the author for
    evaluation.
}
\value{
    This function returns a structure of class \code{ima} composed of the following elements:
    \item{data}{A \code{data.frame} with the actual data.}
    \item{file}{File name.}
    \item{Assay}{Character vector with the assay infromation (platform, S/N and session name).}
    \item{Lot}{Character vector with reagent Lot # and a number of samples.}
    \item{Analytes}{Character vector of analyte names.}
    \item{Units}{Character vector of analyte units.}
    \item{calibrators}{A \code{data.frame} with nominal concentrations of kit calibrators and QCs.}
    \item{Date}{The date the run was processed.}
    \item{Operator}{Name of the analyst that processed the run.}
    \item{Background}{Numeric vector of background values for all analytes.}
}
\note{
    Currently the function can process data files from "Luminex IS" and "xPONENT" software 
    platforms. Tested on run files with few analytes only. Please email problem files to the 
    author.
}
\author{Michal J. Figurski, PhD \email{mfigrs@gmail.com}}
\seealso{\code{\link{sigfit}} and \code{\link{read.kits}}. }
\examples{\dontrun{
    run = read.multiplex("your-path-here","filename.csv")
    run
    plot(run)
    summary(run, analyte=1)
}}
\keyword{ immunoassay }
\keyword{ ima }
\keyword{ read }
\keyword{ csv }
\keyword{ raw data }
\keyword{ Luminex }
\keyword{ xMap }

