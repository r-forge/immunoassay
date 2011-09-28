batch <- function(path, subfolder="", kit.file=NULL, analytes="all",
    run.options     = list(CTS=20, MFI="last", CV=20), 
    model.options   = list(model="L.5", weights="sqrt", refit=0.2, use=1, stvals="adaptive"),
    project.options = list(report="text", template="short", trace=TRUE),
    correct.errors  = NULL) {
        
    # Define internal functions        
    validate <- function(x) {
        # Get list of valid analytes
        v.a      = names(x)[match(paste("pred", attributes(x)$Analytes, sep="."), names(x))]
        analytes = unlist(strsplit(v.a, split=".", fixed=T))[seq(1, length(v.a),1)*2]
        
        # Get lowest standard's name
        standards = unique(as.character(x$SPL[x$Type=="Standard"]))
        lastst    = run.options$MFI
        if (run.options$MFI=="last")  { lastst = standards[length(standards)] }
        if (run.options$MFI=="first") { lastst = standards[1] }
        
        # Validate
        for (i1 in 1:length(analytes)) {
            # Criteria #1 - COUNTS
            x[x[,paste("CTS", analytes[i1], sep=".")] < run.options$CTS, v.a[i1]] = NA
            
            # Criteria #2 - MFI must be greater than that of the lowest standard
            fit       = fits[[i1]]$data
            use       = fit$use[unlist(strsplit(rownames(fit), " ", fixed=T))[seq(1, nrow(fit)*2, 2)]==lastst]
            mfith     = mean(x[x$SPL==lastst, paste("MFI", analytes[i1], sep=".")]*use, na.rm=T)
            if (is.na(mfith)) {
                pos   = match(lastst, standards)
                if (pos>1) { 
                    lastst = standards[pos-1] 
                    use    = fit$use[unlist(strsplit(rownames(fit), " ", fixed=T))[seq(1, nrow(fit)*2, 2)]==lastst]
                    mfith  = mean(x[x$SPL==lastst, paste("MFI", analytes[i1], sep=".")]*use, na.rm=T)
                    }
                }
            if (is.na(mfith)) { mfith = 0 }
            x[x$SPL!=lastst & !is.na(x[,paste("MFI", analytes[i1], sep=".")]) & 
                 x[,paste("MFI", analytes[i1], sep=".")] < mfith, v.a[i1]] = NA
                 
            # Recalculate CVs
            for (i2 in unique(x$SPL)) {
                bz = x[x$SPL %in% i2, v.a[i1]]
                bzz = (sd(bz, na.rm=T)/mean(bz, na.rm=T) * 100)
                if (is.na(bzz)) bzz = 0
                x[x$SPL == i2, paste("pred", analytes[i1], "cv", sep=".")] = bzz
                }  
                 
            # Criteria #3 - CV between replicates must be less than X%
            x[!(x$Type %in% c("Standard","QC")) & !is.na(x[,paste("pred", analytes[i1], "cv", sep=".")]) & 
                (x[,paste("pred", analytes[i1], "cv", sep=".")] > run.options$CV) , v.a[i1]] = NA
            }
        invisible(x)
    }
    
    if (project.options$trace) cat("Checking folders and files ")
    # Get the list of valid files: "files" vector, "path" value, "N" value, "n" value, "l.analytes" vector
    {
    # Check folders first
    if (!file.exists(path=path)) stop("The specified path: \"",path,"\" could not be found.")
    if (nchar(subfolder)>0) 
        if (substr(subfolder,0,1)!="/" & substr(subfolder,0,2)!="\\") subfolder = paste("/", subfolder, sep="")
    if (!file.exists(path=paste(path, subfolder, sep=""))) stop("The specified subfolder: \"",
        subfolder,"\" could not be found within the path: \"", path,"\".")
    ppath  = paste(path, subfolder, sep="")
    files = list.files(path=ppath)
    i1    = grep(".csv", files)
    if (length(i1)>0) files = files[i1]
    else stop("No \".csv\" data files found in the selected folder.")
    ppath  = paste(ppath, "/", sep="")
    
    # Check files if valid Immunoassay runs
    for (i1 in length(files):1) {
        a     = scan(paste(ppath, files[i1], sep=""), what="character", sep=",", nlines=1, quiet=T)
        if (length(grep("xPONENT",a))==0 & length(grep(" IS",a))==0) files = files[-i1]
        if (project.options$trace) cat(".")
        }
    if (project.options$trace) cat(" - done!\n")
    if (length(files)==0) stop("No supported immunoassay data files in the selected folder.")
    else {
        cat("The following valid (supported) immunoassay run files found:\n")
        print(files)
        #cat("Process them? (y/n): ")
        #input = scan(what="character", nmax=1, quiet=T)
        #if (!(input %in% c("y","Y","n","N"))) stop("Undefined choice.")
        #if (tolower(input)=="n") return()
        }

    # Get the kits data: "immunoassay.KITS" vector
    if (!exists("immunoassay.kits")) { 
        immunoassay.kits <<- try(read.kits(path, kit.file))
        if (class(immunoassay.kits) == "try-error") immunoassay.kits <<- try(read.kits(ppath, kit.file))
        if (class(immunoassay.kits) == "try-error") stop("Kit data file:", kit.file, "not found.")
        }

    # Set globals
    N = length(files)
    a = read.multiplex(ppath, files[i1], analytes=analytes)
    l.analytes = attr(a, "Analytes")
    n = length(l.analytes)
    }
    
    if (project.options$trace) cat("Processing starting values list ")
    # Prepare starting value matrix: "COEFS" vector
    {
    coefs.4 = data.frame(matrix(ncol=4, nrow=1, dimnames=list(c("start"), c("a","b","c","d"))))
    coefs.4["start",] = c(a=-100, b=20000, c=100, d=-1)
    coefs.5 = data.frame(matrix(ncol=5, nrow=1, dimnames=list(c("start"), c("a","b","c","d","f"))))
    coefs.5["start",] = c(a=-10,  b=25000, c=200, d=-1.5, f=1)
    coefs   = vector(mode="list", 6); names(coefs) = c("none","123","248","sqrt","1/y","custom")
    coefs   = list(coefs, coefs, coefs, coefs); names(coefs) = c("H.4","H.5","L.4","L.5")
    for (i1 in c(1,3)) { for (i2 in 1:6) { coefs[[i1]][[i2]] = coefs.4 } }
    for (i1 in c(2,4)) { for (i2 in 1:6) { coefs[[i1]][[i2]] = coefs.5 } }
    l.coefs = vector(mode="list", n)
    for (i1 in 1:n) { l.coefs[[i1]] = coefs }
    names(l.coefs) = l.analytes
    if (project.options$trace) cat("- done!\n")
    if (exists("immunoassay.coefs")) {
        
        cat("Data.frame \"immunoassay.coefs\" exists. Use it (u) or overwrite (o)?\n")
        input=scan(what="character", nmax=1, quiet=T)
        if (!(input %in% c("u","U","o","O"))) stop("Undefined choice.")
        if (tolower(input)=="o") immunoassay.coefs <<- l.coefs
        else l.coefs = immunoassay.coefs
        }
    else immunoassay.coefs <<- l.coefs
    }

    if (project.options$trace) cat("Processing model options ")
    # Process options lists
    for (i1 in 1:5) {
        if (project.options$trace) cat(".")
        if (class(model.options[[i1]])=="list") { 
            if (length(model.options[[i1]])==n) {       # if params for all analytes - replicate for all files
                a = vector(mode="list", N)
                for (i2 in 1:N) { a[[i2]] = model.options[[i1]] }
                names(a) = files
                model.options[[i1]] = a
                }
            else { 
                if (length(model.options[[i1]])!=N) stop("List of model parameters must be the 
                    same length as the number of data files.") 
                else {
                    for (i2 in 1:N) { if (length(model.options[[i1]][[i2]])!=n) stop("Parameter \"",
                        names(model.options)[i1], "\" for file \"", files[i2],"\" not defined for all analytes.")
                        }
                    }
                }
            }
        else {
            a  = vector(mode="list", N)
            aa = vector(mode="list", n)                   
            for (i2 in 1:n) { aa[[i2]] = model.options[[i1]] }  # replicate for each analyte
            names(aa) = l.analytes
            for (i2 in 1:N) { a[[i2]]  = aa }                   # replicate for each file
            names(a)  = files
            model.options[[i1]] = a
            }
        }
    if (project.options$trace) cat(" - done!\n")
    # Apply corrections
    if (!is.null(correct.errors)) {
        for (i1 in names(correct.errors)) {
            for (i2 in 1:n) {
                if (correct.errors[[i1]]$name=="use") model.options$use[[i1]][[i2]] = correct.errors[[i1]]$value
                }
            }
        }
    if (exists("immunoassay.options")) {
        if (all(names(immunoassay.options$weights) == files)) {
            cat("Run options already defined in \"immunoassay.options\". Use it (u) or overwrite (o)?\n")
            input = scan(what="character", nmax=1, quiet=T)
            if (!(input %in% c("u","U","o","O"))) stop("Undefined choice.")
            if (tolower(input)=="o") immunoassay.options <<- model.options
            else model.options = immunoassay.options
            }
        else immunoassay.options <<- model.options
        }
    else immunoassay.options <<- model.options
       
    if (project.options$trace) cat("Fitting models: \n")
    if (project.options$report=="text" & file.exists(paste(path, "batch-report.txt", sep="/"))) {
        file.remove(file = paste(path, "batch-report.txt", sep="/")) }
    # For each valid file in the folder
    for (i1 in 1:N) {
        # Read data
        l.run      = read.multiplex(ppath, files[i1], analytes=analytes)
        if (project.options$trace) cat(" + running \"",files[i1],"\" file ")
        if (!all(attr(l.run, "Analytes") %in% l.analytes)) stop("Run files in this folder have different analytes.")
        
        # Fit the data
        fits = vector(mode="list", n)
        for (i2 in 1:n) {
            if (project.options$trace) cat(".")
            zz      = sigfit(l.run, i2, 
                      model   = model.options[["model"]][[i1]][[i2]], 
                      weights = model.options[["weights"]][[i1]][[i2]], 
                      use     = model.options[["use"]][[i1]][[i2]], 
                      refit   = model.options[["refit"]][[i1]][[i2]], 
                      stvals  = model.options[["stvals"]][[i1]][[i2]])
            if (project.options$trace) cat("|")                      
            l.run$Kit   = attr(l.run, "Kit")[1]
            l.run       = predict(zz, l.run, e.fit=T)
            fits[[i2]]  = zz
            }   
        
        # Validate
        if (project.options$trace) cat(".")
        l.run = validate(l.run)
        
        # Create report
        if (project.options$trace) cat(".")
        if (project.options$report=="text") {
            plik = file(description = paste(path, "batch-report.txt", sep="/"), open = "at", blocking=TRUE)
            sink(file = plik, append = TRUE, type = "output")
            if (project.options$template=="full") {
                for (i2 in 1:n) {
                    cat("\n\n# ************************************************************************************** #\n\n")
                    summary(fits[[i2]])
                    }
                cat("\n\n# ************************************************************************************** #\n\n")
                cat("\n\n# ************************************************************************************** #\n\n")
                }
            if (project.options$template=="short") {
                for (i2 in 1:n) {
                    bz  = fits[[i2]]
                    bzm = paste(bz$model[1], bz$model[2], sep=".")
                    bzw = paste(bz$model[3])
                    cat("\nFile:", bz$file, "; Analyte:", bz$analyte[1],"\n")
                    cat("  + model:      ", bzm, "; weights: ", bzw, "\n")
                    cat("  + calibrators:", bz$data$use, "\n")
                    if (any(!is.na(bz$stats))) {
                        cat("  + SSE:        ", bz$stats$SSE[bz$stats$model==bzm & bz$stats$weights==bzw], "\n")
                        cat("  + R-squared:  ", bz$stats$r.sq[bz$stats$model==bzm & bz$stats$weights==bzw], "\n")
                        }
                    }
                }
            sink()
            close(plik)
            }
        if (project.options$report=="latex") {
            immunoassay.environment <<- environment()
            if (!file.exists(path=paste(path, "Tex", sep="/"))) dir.create(paste(path, "/Tex", sep=""))
            if (!file.exists(path=paste(path, "Tex/pdf", sep="/"))) dir.create(paste(path, "/Tex/pdf", sep=""))
            Sweave(file=paste(path, "/", project.options$template, ".run-report.r", sep=""), output=paste(path, "/Tex/", 
                attr(l.run, "file"),".tex", sep=""), quiet=T)
            if (project.options$trace) cat("*")
            }
        
        # Collect results
        if (project.options$trace) cat(".")        
        if (exists("results")) { results = rbind(results, as.data.frame(l.run)) }
        else results = as.data.frame(l.run)
        
        # The end
        if (project.options$trace) cat("\n")
        }
    if (project.options$report=="latex") {
        if (project.options$trace) cat("* running final report\n")
        if (file.exists(paste(path, "/", project.options$template, ".project-report.r", sep=""))) {
            Sweave(file=paste(path, "/", project.options$template, ".project-report.r", sep=""), 
                output=paste(path, "/Tex/", project.options$template, ".project-report.tex", sep=""), 
                quiet=T)
            }
        else warning("Project report template file not found.")
        }

    return(results)   
}

