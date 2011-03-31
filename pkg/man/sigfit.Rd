\name{sigfit}
\alias{sigfit}
\title{Fit the sigmoidal model to raw immunoassay run data.}
\description{
    This function fits either a 4-PL or 5-PL sigmoidal model (classic logistic model and the 
    Hill's form) to the data from a immunoassay run - nominally an object of class 
    \code{"ima"}. The user can select the model type, several types of weighting, 
    which calibrators to use and even the starting values for the \code{nls} fit. This 
    function can also automatically select best fit from user selected list, and it can 
    remove erroneous calibration points (only one replicate) automatically to improve the fit.
}
\usage{
    sigfit(x, analyte = 1, model = "L.5", weights = "sqrt", 
        refit = 0.2, use = 1, stvals = "adaptive", sledz=NULL)
}
\arguments{
  \item{x}{An object of class \code{"ima"}.} 
  \item{analyte}{
    Numeric. Selects the analyte for which the fit is desired - from the analytes
    available in the data
    }
  \item{model}{
    Character vector, can be any combination of the following values: "L.4", "L.5", 
    "H.4" and "H.5", or alternatively it can have the value \code{"auto"} which is 
    equivalent to providing all four model types. "L" represent classic Logistic equation models, 
    while "H" represent Log-Logistic (Hill form) models. "4" and "5" after the dot represent "4-PL" 
    and "5-PL" models, respectively. Defaults to "L.5", which we found to be the most robust model.
    }
  \item{weights}{
    Character or Numeric vector. If \code{character} - it can be the name or combination 
    of names of any of the five built-in weighting types, or the value of \code{"auto"} 
    which is equivalent to providing a vector of all weighting type names. The weighting
    types are as follows: "1/y" - standard 1/y weighting; "sqrt" - square root of 1/y
    weights; "248" - calibrators are assigned weights of consecutive powers of 2 starting 
    from the lowest calibrator; "123" - similarly to "248", calibrators are assigned 
    weights of consecutive natural numbers; "none" - no weighting. \cr
    If \code{numeric} - it must be a vector of weights of the same length as the number of
    calibrator replicates. \cr
    Defaults to "sqrt", which according to our experience is the most robust weighting type.
    }
  \item{refit}{
    Numeric. If automatic selection of calibrator replicates is desired, this is the threshold
    for the in-accuracy error of the standard. If the post-hoc estimated value for any of the 
    replicates has error above the given threshold, the function will attempt to re-fit the 
    model without this replicate and check if this improves the fit. This can be suppressed 
    by providing value of \code{0} or \code{NA}. Defaults to 0.2 (20\% in-accuracy). \cr
    Calibrators removed by this criteria show as "x" on the plot made with the 
    \code{\link{plot.sigfit}} function.
    }
  \item{use}{
    A numeric vector of length one, or the exact same length as the number of calibrator 
    replicates. It is to be composed of "0"s and "1"s and can be used to arbitrarily select 
    calibrator replicates to use for fitting. Use the value of \code{1} to turn the replicate 
    "on" or the value of \code{0} to turn it "off". Calibrators turned "off" manually show 
    as small solid dots on the plot made with the \code{\link{plot.sigfit}} function.\cr
    Defaults to 1, which means to use all replicates.
  }
  \item{stvals}{
    Character or a list. It can have the following values: NULL - tells the function to use
    the default starting values for the \code{nls} fit; "adaptive" (\code{class "character"})
    - tells the function to use adaptive starting values - useful when fitting more than 
    just a few runs, the function will use medians of the coefficients of already fitted models
    for each analyte to derive starting values for subsequent fits. This can sometimes 
    stabilize the behaviour of the fitting function in larger projects, and this is the 
    default.\cr
    If a list, it must be a named list of starting values of all model parameters with their
    values: \code{"list(a=, b=, c=, d=)"} for 4-PL models and \code{"list(a=, b=, c=, d=, f=)"}
    for 5-PL models. If using the list for this parameter, currently the function will 
    allow only one model type (either 4-PL or 5-PL), but multiple weighting types are 
    still possible.
    }
  \item{sledz}{
    Logical. This turns on internal debugging for use by the developer only. Defaults to FALSE.
    }
}
\details{    
    The \code{sigfit} function is the actual work horse of the whole \code{immunoassay} library. 
    It does the fitting of sigmoidal models by using the \code{nls} function. Four models
    are hard-coded in it to choose from: two 4-PL and two 5-PL models, each can be either 
    in the form of standard sigmoidal equation (logistic function), or in the logarithm (Hill's)
    form.
    
    The user can either select the model and all its parameters explicitly, or can provide
    vectors of parameters for the function to choose from. The function will fit the models
    for all combinations of parameters and will select the best-fit model based on the 
    criteria of residual standard error and R-squared. The user can also set the model and
    weighting parameters to \code{"auto"} and leave the decision on selecting the model 
    entirely in the hands of the algorithm in this function. 
    
    The function can attempt to recognize and remove single replicates of standards that 
    cross the threshold of in-accuracy (\code{refit} parameter). Of note: once the threshold is 
    crossed by at least one replicate, the function will go through all replicates in the 
    sequence of decreasing error (in-accuracy) and attempt to re-fit the model without 
    them. This may sometimes result in more than one replicate removed (but always only
    one replicate per standard), even if only one replicate was above the threshold in the
    beginning. The criteria to "keep the replicate out" is currently that the overall 
    error of the old model must be more than \code{1+refit} times the error for the 
    re-fit model. Additionally, the QC criteria are checked: whether the QC predictions
    for the re-fit model are within the limits - if they were within those limits for the
    old model. If not, the removed standard is re-introduced - this is to prevent situations
    where automatic removing of calibrators would improve the fit but mess up the QCs.
    
    In situations where the user is certain that some standards must be excluded from fitting,
    it can be done using the \code{use} parameter - simply put "0"s for the replicates to
    be removed. 
    
    If the model convergence fails with the default \code{nls} starting values, that are 
    hard-coded, users are encouraged to experiment with the \code{"stvals"}. The starting 
    values can be provided directly as a list, or the option \code{"adaptive"} can be used
    as well. For the \code{"adaptive"} starting values to work, a global list named 
    \code{"immunoassay.coefs"} must be present and appropriately formatted. This list can be
    created manually, though the function \code{\link{batch}} creates it automatically if
    it does not exist.   
}
\value{
    Function \code{sigfit} returns an object of class \code{"sigfit"}, that is a list 
    composed of the following items:
    \item{fit }{The actual \code{nls} fit object.}
    \item{data }{The \code{data.frame} with calibrators data, that was used for the fit.}
    \item{qcs }{Similarly, a \code{data.frame} with Quality Control samples data.}
    \item{model }{Character vector of the final model information, describing the model type 
            and weighting type.}
    \item{analyte }{Character vector of analyte information: first, analyte names, then units.}
    \item{file }{Character, the name of the raw data file, for which the fit was done.}
    \item{stats }{A matrix of summary statistics of the models, if more than one model or
            weighting type was provided as an input.}
}
\author{
    Michal J. Figurski, PhD \email{mfigrs@gmail.com} of the Biomarker Research Laboratory, 
    University of Pennsylvania, Philadelphia, PA.
}
\seealso{\code{\link{immunoassay} package}, \code{\link{batch} function} and
    \code{\link{plot.sigfit}}.}
\examples{\dontrun{
    run = read.multiplex("your-path-here")
    fit = sigfit(run)
    fit
    summary(fit)
    plot(fit)
}}
\keyword{ immunoassay }
\keyword{ sigmoidal }
\keyword{ nls }
\keyword{ sigfit }
