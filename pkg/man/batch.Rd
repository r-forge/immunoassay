\name{batch}
\alias{batch}
\title{Batch-processing function for immunoassay data.}
\description{
    This is a wrapper function making use of other lower-level functions in this package. It looks
    for immunoassay data files in a given folder, processes these files and accumulates results into a 
    single data frame. Optionally, it can create a project report - either in a simple text form,
    or in LaTeX format. 
}
\usage{
batch(path, subfolder = "", kit.file = NULL, analytes = "all", 
    run.options     = list(CTS = 20, MFI = "last", CV = 20), 
    model.options   = list(model = "L.5", weights = "sqrt", 
                      refit = 0.2, use = 1, stvals = "adaptive"), 
    project.options = list(report = "text", template = "short", 
                      trace = TRUE), 
    correct.errors  = NULL)
}
\arguments{
  \item{path}{Character string. Path to the main project folder. A csv file with kit data 
    should be located in that folder. Also, all subfolders are created in this folder.}
  \item{subfolder}{Character string. If the csv files with immunoassay data are located in the 
    subfolder of the main folder, provide its name here. Otherwise, leave blank.}
  \item{kit.file}{Character string. File name of the csv file containing kit information. 
    If a \code{data.frame} named \code{"immunoassay.kits"} exists, it will be used, and this 
    parameter will be ignored.}
  \item{analytes}{Character or numeric. This parameter is passed to \code{\link{read.multiplex}}
    function. See more information there.}
  \item{run.options}{List of 3 elements that provide validation criteria: \code{CTS} - 
    numeric, is the threshold for COUNTS; \code{MFI} - character, either "last" or "first",
    that is the last of first calibrator becomes the criteria for lowest MFI; and \code{CV} 
    - numeric, percent CV that is threshold for between-replicate variability.}
  \item{model.options}{A list of 5 elements. Element \code{model} - a character or list of 
    character values of length N, where "N" is the number of data files in the folder. Element 
    \code{weights} - a character or list of character values of length N. Element \code{refit} - 
    a numeric value or list of numeric values of length N. Element \code{use} - a numeric 
    value or list of numeric values of length N. Element \code{stvals} - a character or 
    list of character values of length N. For more information on these parameters, see help 
    for \code{\link{sigfit}}}
  \item{project.options}{A list of 3 elements. Element \code{report} - character, either "text"
    or "latex" options are currently available, defaults to "text". Element \code{template} 
    - character. For \code{report="text"} option two templates are built-in: \code{"short"}, that 
    produces quick summary of each analyte for each plate, and \code{"full"} - an extended template 
    with more information. Defaults to \code{"short"}. For \code{report="latex"} option, this is the 
    name of a LaTeX template, located in the main project folder. For more information, see 
    details below. Element \code{trace} - logical, debugging option. If true, trace information 
    will be displayed when the function is processing data files. Defaults to \code{FALSE}.}
  \item{correct.errors}{A list - an alternative means to introduce corrections to the fit
    parameters for a small number of items. In this list, each named element is a list. The 
    named element's name must be the name of the datafile that is to be corrected. There 
    must be two elements in each named sub-list: \code{name}, that names the parameter to 
    be corrected (currently only \code{name="use"} is implemented), and the second element 
    \code{value} is the new value of the named parameter. Defaults to \code{NULL}.}
}
\details{The \code{batch} function is the core function of this library. It implements 
    automation in processing of entire folders of immunoassay data files. The function 
    first looks for immunoassay files (currently limited to multiplex data files) under 
    the provided path - it does so intelligently, so it
    can distinguish immunoassay run files from other .csv files. Then it creates a list of names 
    of these files and processes through this list according to provided options. 
    
    For each file from the list, this function loads it using \code{\link{read.multiplex}} and 
    fits the sigmoidal model using \code{\link{sigfit}} function with the set of parameters
    given in \code{model.options}. Next, it validates the fitted data using the criteria 
    provided in \code{run.options} parameter, and finally, it creates a report as set in 
    the \code{project.options}. 
    
    There are two major advantages of using this function, instead of manually processing the 
    data: it can save substantial amount of time in processing of large number of data files,
    while preserving considerable flexibility in setting fit parameters and applying corrections.
    The second major advantage is the ability to create reports. Text reports are useful for
    obtaining general insight into the data and fitting process, but this function can really
    shine in connection with \code{Sweave} LaTeX report templates, that can be elaborate 
    programs on themselves.
    
    In order to use LaTeX template option, the user must create a minimum of two template 
    files. These must be named: "my-project-name.run-report.r" - a template for each plate, and 
    "my-project-name.project-report.r" - a project-wide template incorporating (or not) individual
    plate reports. The \code{"my-project-name"} part of the template name must be provided
    as the \code{template=} parameter, when \code{report="latex"} option is used. These templates
    are "Sweaved" in the process of running this function. 
    
    There is one known problem associated with use of the \code{Sweave} function - it's that the
    \code{Sweave} works in global environment and doesn't "see" the environment within the
    \code{batch} function it has been called from. To overcome this, the \code{batch} function 
    sets a global variable - a pointer to its environment, named \code{"immunoassay.environment"}.
    Objects and data from within \code{batch} function can be then accessed within the 
    \code{Sweave} templates using \code{immunoassay.environment$} pre-fix. 
    
    For LaTeX template programmers, the list of the accessible objects within \code{immunoassay.environment}
    is as follows: \code{validate} - validation function; \code{immunoassay.coefs} (global) - list of all 
    fit coefficients; \code{immunoassay.kits} (global) - a \code{data.frame} with kits information;
    \code{ppath} - full path to the data files; \code{files} - list of valid 
    immunoassay files; \code{N} - number of files (\code{length(files)}); \code{l.analytes} - 
    vector of analyte names; \code{n} - number of analytes (\code{length(l.analytes)});
    \code{l.run} - a \code{data.frame} of class \code{ima}, containing the validated run 
    data; and \code{fits} - a list of fitted \code{nls} models, of length "n".
    After all the plates are processed, additional object becomes available: \code{results} -
    a \code{data.frame} being a collection of all results from all plates, for use in the
    project-wide report. 
}
\value{
  This function returns a simple data frame that is the collection of results from all processed
  files. It does not contain any fit information - if fit information are required, some form
  of report must be called.\cr
  In addition to the returned \code{data.frame} this function by default creates three global 
  objects: the \code{immunoassay.coefs} list - a list of model coefficients for each fitted
  run file, and each model type; the \code{immunoassay.options} - a list of model parameters 
  for all run files; and \code{immunoassay.environment} - a pointer to the environment of the
  \code{batch} function for use by the \code{Sweave} report templates.
}  
\note{
    Examples of LaTeX report templates from an example \code{\link{project}} are located in 
    the "templates" directory in the main package folder. See help for the example project for 
    more information.
}
\author{Michal J. Figurski, PhD \email{mfigrs@gmail.com} of the Biomarker Research Laboratory, 
University of Pennsylvania, Philadelphia, PA.}
\seealso{\code{\link{immunoassay}}  package, example \code{\link{project}}, \code{\link{read.multiplex}}, 
    \code{\link{read.kits}} and \code{\link{sigfit}}.}
\examples{\dontrun{

a <- batch("path-to-my-project", subfolder="Results", 
     kit.file="kit-file-name.csv",  analytes=1:2, 
     model.options   = list(model="L.5", weights=c("sqrt","1/y"), refit=0.2, 
                       use=1, stvals="adaptive"),
     project.options = list(report="latex", template="my-project-name", 
                       trace=FALSE),
     correct.errors  = list(`Plasma Pl072.csv`=list(name="use", 
                       value=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1))))

}}
\keyword{ immunoassay }
\keyword{ sigmoidal }
\keyword{ batch processing }
\keyword{ report }
