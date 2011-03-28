sigfit <- function(x, analyte=1, model="L.5", weights="sqrt", refit=0.2, use=1, stvals="adaptive", sledz=NULL) {
    # Error check & variable set
    if (class(x)[1]!="ima") stop("\"x\" must be of class \"ima\".")
    if (is.na(attr(x, "Analytes")[analyte])) stop("No such analyte in the data.")
    if (!(all(toupper(model) %in% c("H.4","H.5","L.4","L.5","AUTO")))) 
        stop("Unsupported model type(s): ", model[!(toupper(model) %in% c("H.4","H.5","L.4","L.5","AUTO"))])
    if (length(refit)!=1 | !(class(refit) %in% c("logical","numeric"))) stop("Refit must either be \"NA\" or a number in range 0 to 1.")
    if (!(is.na(refit))) if (refit<0 | refit>1) stop("Refit must either be \"NA\" or a number in range 0 to 1.")
    
    # Disambiguation
    if (any(model=="auto") | length(model)>1 | any(weights=="auto") | 
       (all(class(weights)=="character") & length(weights)>1) | !is.na(refit)) 
       sigfit.auto(x, analyte, model, weights, refit, use, stvals, sledz)
    else sigfit.default(x, analyte, model, weights, use, stvals)
}

sigfit.default <- function(x, analyte, model, weights, use, stvals) {
    # Error check & variable set
    if (class(x)[1]!="ima") stop("\"x\" must be of class \"ima\".")
    z.analyte = attr(x, "Analytes")[analyte]
    if (is.na(z.analyte)) stop("No such analyte in the data.")
    type  = unlist(strsplit(model, split=".", fixed=T))[2]
    if (!(type %in% c("4","5"))) stop("Equation of type \"", type, "\" not supported.")
    model  = toupper(unlist(strsplit(model, split=".", fixed=T))[1])
    if (!(substr(model,0,1) %in% c("H","L"))) stop("Model equation \"", model, "\" not supported.")

    # Save QCs
    qcs = x[x$Type=="QC", ]
    if (nrow(qcs)>0) {
        rownames(qcs) = paste(qcs$SPL, " (", unlist(strsplit(as.character(qcs$Loc), 
            split = "(", fixed = T))[seq(2, nrow(qcs), 2)], sep = "")
        qcs         = as.data.frame(qcs[,paste(c("MFI", "conc"), rep(z.analyte,2), sep=".")])
        names(qcs)  = c("MFI","value")
        }
    
    # Get calibrators
    x = x[x$Type=="Standard",]
    cals = x[tolower(x$SPL)!="background", paste(c("MFI","conc"), rep(z.analyte,2), sep=".")]
    rownames(cals) = paste(x$SPL, " (", unlist(strsplit(as.character(x$Loc), 
    split = "(", fixed = T))[seq(2, nrow(x), 2)], sep = "")
    cals$weights = 1
    if (length(use)!=1 & length(use)!=nrow(cals)) stop("\"use\" vector must be of length 1 or length equal to number of calibrators.")
    cals$use     = use

    # Get starting values
    if (any(stvals=="adaptive")) {
        if (exists("immunoassay.coefs")) {
            st = immunoassay.coefs[[analyte]][[model]][[weights]]
            if (!is.null(st)) {
                st = na.omit(st) 
                stvals  = if (type=="4") {
                    list(a=median(st[,1])-10, b=median(st[,2])*0.9, c=median(st[,3])*0.9, d=median(st[,4])*0.9) }
                    else {
                    list(a=median(st[,1])-10, b=median(st[,2])*0.9, c=median(st[,3])*0.9, d=median(st[,4])*0.9, 
                        f=median(st[,5])*0.9) }
                    }
                else stvals = immunoassay.coefs[[analyte]][[model]][[weights]]["start",]
                }
            else stvals = NULL
            }
    if (is.null(stvals)) { 
        stvals = list(c(a=-100, b=20000, c=100, d=-1), c(a=-100, b=20000, c=100, d=-1,
            f=1))[[ifelse(type=="4",1,2)]] 
        if (toupper(substr(model,0,1))=="H") stvals["c"] = log(stvals["c"])
        }
    if (class(stvals)!="numeric") stop("Check your starting values.")

    # Set weights
    for (i1 in unique(cals[,2])) cals$aveMFI[cals[,2]==i1] = mean(cals[cals[,2]==i1,1], na.rm=T)
    if (length(weights)==0) { weights = "none" }
    if (length(weights)==1) if (is.na(weights)) { weights = "none" }
    if (length(weights)==1 & class(weights)=="character") {
        if (!(weights %in% c("1/y","sqrt","248","123","none"))) stop("Unknown weighting type.")
        if (weights=="1/y")  cals$weights = 1/(cals$aveMFI^2)/(max(1/(cals$aveMFI^2), na.rm=T))
        if (weights=="sqrt") cals$weights = sqrt(1/(cals$aveMFI^2)/(max(1/(cals$aveMFI^2), na.rm=T)))
        if (weights=="248")  cals$weights = rev(1/(2^rep(0:(floor(nrow(cals)/2)-1), each=2)))
        if (weights=="123")  cals$weights = rev(1/(rep(1:(floor(nrow(cals)/2)), each=2)))
        if (weights=="none") cals$weights = rep(1, nrow(cals))
        }
    else {
        if (length(weights)!=nrow(cals)) stop("Length of \"weights\" must be the same as the number of calibrators (", nrow(cals),").")
        if (!(class(weights) %in% c("numeric","integer"))) stop("Unknown weighting type.")
        if (any(weights>1)) weights = weights / max(weights, na.rm=T)
        cals$weights = weights
        }

    # Apply "use" to weights:
    cals$weights[(!is.na(cals[,1]) & cals[,1]<=0) | (!is.na(cals$use) & cals$use<=0)] = 0
    cals$weights[is.na(cals[,1]) | is.na(cals$use)] = 0
    names(cals)[1:2] = c("MFI","value")
    cs      = na.omit(cals)

    # Do fitting
    if (model=="H") { 
        if (type=="4") z.fit = nls(MFI ~ a + b/(1+10^((c-log(value))*d)), 
            data    = cs, 
            start   = stvals, 
            control = list(maxiter=2000, warnOnly = F) , algorithm="port",
            weights = cs$weights)
        if (type=="5") z.fit = nls(MFI ~ a + b/(1+10^((c-log(value))*d))^f, 
            data    = cs, 
            start   = stvals, 
            control = list(maxiter=2000, warnOnly = F) , algorithm="port",
            weights = cs$weights)
        }
    if (model=="L") {
        if (type=="4") z.fit = nls(MFI ~ a + b/((1+(value/c)^d)), data=cs, 
            start   = stvals, 
            weights = cs$weights, 
            control = list(maxiter=1000, warnOnly = F) , algorithm="port")
        if (type=="5") z.fit = nls(MFI ~ a + b/((1+(value/c)^d)^f), data=cs, 
            start   = stvals, 
            weights = cs$weights, 
            control = list(maxiter=1000, warnOnly = F) , algorithm="port")
        }

    # Finish
    return(structure(list(
        fit     = z.fit, 
        data    = as.data.frame(cals)[1:4],
        qcs     = qcs,
        model   = c(equation = model, type=type, weighting=ifelse(class(weights)=="character",weights,"custom")), 
        analyte = c(analyte=z.analyte, unit=attr(x, "Units")[analyte]),
        file    = attr(x, "file"),
        stats   = NA), 
        class   = "sigfit"))
}

sigfit.auto <- function(x, analyte, model, weights, refit, use, stvals, sledz) {
    # Set model lists and check errors
    if (class(x)[1]!="ima") stop("\"x\" must be of class \"ima\".")
    if (is.na(attr(x, "Analytes")[analyte])) stop("Selected analyte does not exist in the data.")
    if (class(model)=="character") {
        if (any(tolower(model)=="auto")) m.list=c("H.4","H.5","L.4","L.5")
        else {
            if (all(toupper(model) %in% c("H.4","H.5","L.4","L.5"))) m.list=toupper(model)
            else stop("Unsupported model type(s): ", model[!(toupper(model) %in% c("H.4","H.5","L.4","L.5"))])
            }
        }
    else stop("Unsupported model type(s): ", model[!(toupper(model) %in% c("H.4","H.5","L.4","L.5"))])

    # Set weight list & check for errors
    if (class(weights)=="character") {
        if (any(tolower(weights)=="auto")) w.list=c("none","123","248","sqrt","1/y")
        else {
            if (all(tolower(weights) %in% c("none","123","248","sqrt","1/y"))) w.list=tolower(weights)
            else stop("Unsupported weighting type(s): ", weights[!(tolower(weights) %in% c("none","123","248","sqrt","1/y"))])
            }
        }
    else { w.list = "custom" }
    
    # Get all fits
    ssr  = matrix(ncol=11, nrow=0)
    fits = vector(mode="list", length=length(m.list)*length(w.list)); counter=1
    for (i2 in m.list) {
        for (i3 in w.list) {
            if (!is.null(sledz)) { cat("\n", sledz, "*******", i2, i3, "*******")
            cat("\n", sledz, "* Podejscie 1") }  # Po prostu wpasuj funkcje
            
            if (i3!="custom") z.fit = try(sigfit.default(x, analyte, model=i2, weights=i3, use, stvals), silent=T)
            else z.fit = try(sigfit.default(x, analyte, model=i2, weights=weights, use, stvals), silent=T)
            if (class(z.fit)=="try-error") {
                if (!is.null(sledz)) { cat("\n", sledz, "* Podejscie 2") } # Moze uzyj standardowych wartosci startowych
                if (length(grep("nls", z.fit))==0 & length(grep("numericDeriv", z.fit))==0) stop(geterrmessage())
                if (i3!="custom") z.fit = try(sigfit.default(x, analyte, model=i2, weights=i3, stvals=NULL, ...), silent=T)
                else z.fit = try(sigfit.default(x, analyte, model=i2, weights=weights, stvals=NULL, ...), silent=T)
                if (class(z.fit)=="try-error") {
                    if (!is.null(sledz)) { cat("\n", sledz, "* Podejscie 3") } # Wywal 1szy lub ostatni kalibrator
                    if (length(use)==1) uzyj = rep(1,nrow(x[x$Type=="Standard",]))
                    else uzyj = use
                    temp = x[x$Type=="Standard", c("SPL", paste("MFI", attr(x, "Analytes")[analyte], sep="."))]
                    cals = as.character(temp$SPL)
                    if (any(temp[cals %in% cals[length(cals)], 2] < 10 * attr(x, "Background")[analyte])) {
                        uzyj[cals %in% cals[length(cals)]] = NA }
                    else { uzyj[cals %in% unique(cals)[1]] = NA }
                    if (i3!="custom") z.fit = try(sigfit.default(x, analyte, model=i2, weights=i3, stvals=NULL, use=uzyj), silent=T)
                    else z.fit = try(sigfit.default(x, analyte, model=i2, weights=weights, stvals=NULL, use=uzyj), silent=T)
                    if (class(z.fit)=="try-error") {
                        if (!is.null(sledz)) { cat("\n", sledz, "* Nie udalo sie.") }
                        ssr   = rbind(ssr, c(analyte, i2, i3, rep(NA,8)))
                        counter = counter+1
                        next
                        }
                    else { if (!is.null(sledz)) { t1 = check(z.fit); cat("\n", sledz, "* SSE:", t1$SSE, "; sigma:", t1$sigma, "; R-kwadrat:",t1$r.squared) } }
                    }
                else { if (!is.null(sledz)) { t1 = check(z.fit); cat("\n", sledz, "* SSE:", t1$SSE, "; sigma:", t1$sigma, "; R-kwadrat:",t1$r.squared) } }
                }
            else { if (!is.null(sledz)) { t1 = check(z.fit); cat("\n", sledz, "* SSE:", t1$SSE, "; sigma:", t1$sigma, "; R-kwadrat:",t1$r.squared) } }
            pred = try(predict(z.fit, e.fit=T), silent=T)
            if (!is.null(sledz)) { cat("\n", sledz, "***********************") }
            if (class(pred)=="try-error") {
                ssr    = rbind(ssr, c(analyte, i2, i3, rep(NA,8)))
                counter  = counter+1
                warning(paste("Fit", i2, i3, "for:", analyte, "unreliable"))
                next
                }
            z.cals = levels(factor(x$SPL[x$Type=="Standard"]))
            ssr    = rbind(ssr, c(analyte, i2, i3, unlist(check(z.fit))))

            # Global store
            if (exists("immunoassay.coefs")) immunoassay.coefs[[analyte]][[i2]][[i3]][attr(x, 
                "file"), ] <<- summary(z.fit$fit)$coef[,1]
            
            # Keep model object
            fits[[counter]] = z.fit
            counter         = counter+1
            }
        }
    
    # Select the best fit (sigma / R-squared criteria)
    ssr   = data.frame(ssr)
    names(ssr) = c("analyte","model","weights","st.err.median","st.err.mean","qc.err.median",
        "qc.err.mean","SSE","sigma","Syx","r.squared")
    for (i4 in 4:ncol(ssr)) {
        ssr[,i4] = as.numeric(as.character(ssr[,i4]))
        if (is.na(refit)) { ssr[!is.na(ssr[,i4]) & ssr[,i4]<0,i4] = NA }
        }    
    ssr$criteria   = (ssr$sigma/ssr$r.squared)
    if (length(m.list)==1 & length(w.list)==1) {
        z.fit=fits[[1]]
        if (is.null(z.fit)) stop("Fit for model \"", model, "\" with \"", 
            ifelse(class(weights)=="character",weights,"custom"), "\" weighting for ", 
            attr(x, "Analytes")[analyte]," failed.")
        by.all=1
        }
    else {
        repeat {
            by.all    = suppressWarnings(match(min(abs(ssr$criteria), na.rm=T), abs(ssr$criteria)))
            z.fit  = fits[[by.all]]
            if (is.null(z.fit)) stop("No reliable fit could be automatically obtained from the provided list of models.")
            if (any(is.na(predict(z.fit))) & is.na(refit)) { ssr$criteria[by.all] = NA  } # Eliminate models that make NA predictions in any calibrator
            else break
            }
        }
    #z.fit       = fits[[by.all]]
    if (!is.null(sledz)) { t1 = z.fit$model; cat("\n", sledz, "=> Koncowy:", paste(t1[1],t1[2], sep="."), "; wagi:",t1[3],"\n") }
    z.fit$stats = ssr[,-c(1,12)]
    
    # If necessary, remove calibrators and re-fit
    if (!is.na(refit) & refit>0) {
        if (!is.null(sledz)) { cat("\n", sledz, "+ Refit") }
        if (length(use)==1) uzyj = rep(1,nrow(x[x$Type=="Standard",]))
        else uzyj = use
        use.old = uzyj 
        usen    = as.character(x[x$Type=="Standard","SPL"])
        pred    = predict(z.fit, e.fit=T)
        if (any(abs(pred$error)>refit*100, na.rm=T)) {
            # Calculate rank of error
            pred$rank[order(abs(pred$error), decreasing=T)] = seq(1,nrow(pred),1)  
            
            # Go through all ranks, starting from the biggest error
            for (i5 in 1:(max(na.omit(pred$rank))-2)) {
                n  = match(i5, pred$rank)                # Find which calibrator has the rank
                if (!is.null(sledz)) { cat("\n", sledz, "+ Rank:", i5) }
                # If one replicate is already NA or 0, skip:
                if (any(is.na(uzyj[usen %in% usen[n]])) | any(uzyj[usen %in% usen[n]]==0)) { next }
                uzyj[n] = NA                                        # Make the calibrator NA
                if (!is.null(sledz)) { cat(" uzyj: ", uzyj) 
                n.fit  = try(sigfit(x=x, analyte=analyte, model=model, weights=weights, 
                    refit=NA, use=uzyj, stvals="adaptive", sledz=paste(sledz,"   ")), 
                    silent=T)                                       # Refit with new calibrators
                    }
                else {
                n.fit  = try(sigfit(x=x, analyte=analyte, model=model, weights=weights, 
                    refit=NA, use=uzyj, stvals="adaptive", sledz=NULL), silent=T)
                    }
                if (class(n.fit)=="try-error") { # If refit failed
                    uzyj[n] = use.old[n]                            # Restore calibrator
                    if (!is.null(sledz)) { cat("\n", sledz, "- Re-fit sie nie powiódl.\n") }
                    next }# Skip the rest
                val.old = check(z.fit)                  # Evaluate the fit
                val.new = check(n.fit)

                if (!is.null(sledz)) { cat("\n", sledz, "- Re-fit udany. Stare sigma:", 
                    val.old$sigma, "; nowe sigma*", 1+refit, ": ", val.new$sigma*(1+refit)) }

                # Now check the refit
                # If error over 500%, don't even check - remove it!
                if (!is.na(pred$error[n])) if (abs(pred$error[n])>500) { 
                    if (!is.null(sledz)) { cat("\n", sledz, "- Blad wiekszy niz 500%") }
                    z.fit=n.fit; next }   
                # Check sigma first 
                if (exists("t1")) rm(t1)
                if (val.old$sigma > (1+refit) * val.new$sigma) {
                    if (!is.null(sledz)) { cat("\n", sledz, "- Nowe sigma*", 1+refit, " mniejsze od starego") }
                    t1 = as.data.frame(predict(n.fit, newdata=x[x$Type=="QC",]))
                    t2 = attributes(x)$qc.ranges
                    t3 = as.character(unique(x$SPL[x$Type=="QC"]))
                    t0 = TRUE
                    for (i6 in 1:length(t3)) {
                        t4  = t1[t1$SPL==t3[i6], paste("pred",attr(x,"Analytes")[analyte], sep=".")]
                        t5l = paste("Con", LETTERS[i6], ".lo", sep="")
                        t5h = paste("Con", LETTERS[i6], ".hi", sep="")
                        if (!is.null(sledz)) { cat("\n", sledz, "-", substr(t5l,0, nchar(t5l)-3), 
                            "od", t2[t2$Analyte==attr(x,"Analytes")[analyte], t5l], 
                            "do", t2[t2$Analyte==attr(x,"Analytes")[analyte], t5h], "wartosci:", t4 ) }
                        if (any(t4 < t2[t2$Analyte==attr(x,"Analytes")[analyte], t5l]) | 
                            any(t4 > t2[t2$Analyte==attr(x,"Analytes")[analyte], t5h])) t0=FALSE
                        }
                    if (!is.null(sledz)) { cat("\n", sledz, "- Kryterium zakresu QC:", ifelse(t0," przeszlo", " padlo")) }
                    if (!t0) {
                        t4 = t1[, paste(c("conc","pred"), rep(attr(x,"Analytes")[analyte],2), sep=".")]
                        names(t4) = c("val","pred")
                        if (all(abs(((t4$pred-t4$value)/t4$value))<refit)) t0=TRUE
                        if (!is.null(sledz)) { cat("\n", sledz, "- Poprawka kryterium zakresu QC:", ifelse(t0," przeszlo", " padlo")) }
                        }
                    if (t0) { z.fit = n.fit; next }                    
                    }
                # Then check if QCs are improved
                if (!is.na(val.old$QC.error["median"]) & !is.na(val.old$QC.error["median"])) {
                    if (val.old$QC.error["median"] > 1.1 * val.new$QC.error["median"]) {
                        if (val.old$St.error["median"] < 1.1 * val.new$St.error["median"]) { break }
                        z.fit = n.fit
                        if (!is.null(sledz)) { cat("\n", sledz, "- Kryterium mediany QC zastosowane") }
                        next }
                    else {
                        uzyj[n] = use.old[n]
                        if (!is.null(sledz)) { cat("\n", sledz, "- Kryterium mediany QC nie przeszlo\n") }
                        next 
                        }
                    }
                }
            }
        else if (!is.null(sledz)) { cat(" - ... refit criteria not met.\n") }

        #if (length(model)==1 & length(weights)==1) {
        #    ssrx = matrix(ncol=10, nrow=0)
        #    ssrx = as.data.frame(rbind(ssrx, c(model=model, weights=weights, unlist(check(z.fit)))))
        #    for (i6 in 3:ncol(ssrx)) ssrx[,i6] = as.numeric(as.character(ssrx[,i6]))
        #    z.fit$stats = ssrx
        #    }
        }

    # Return fit data
    invisible(z.fit)
}

