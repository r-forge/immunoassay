plot.sigfit <- function(x, table=TRUE, type="fit", norm="weighted", ...) {
    require(plotrix)
    data = x$data
    analyte = x$analyte["analyte"]
    units   = x$analyte["unit"]
    if (is.na(units)) units = ""
    symbol = data$use+20
    symbol[is.na(symbol)] = 4
    if (type=="fit") {
        plot(MFI~value, data=data, pch=symbol, xlab=paste(analyte, "concentration", units), 
            ylab=paste(analyte,"MFI"),  ...)
        newdata   = data.frame(value=seq(min(data$value, na.omit=TRUE)*0.1, max(data$value, na.omit=TRUE), length=50))
        newdata$y = predict(x$fit, newdata=newdata)
        lines(newdata, col="#ff121299", lwd=2)
        if (table) addtable2plot(0,max(data$MFI)*0.85, data.frame(coefficient=formatC(coefficients(x$fit), digits=9, format="fg")),
            bty="n", cex=0.6, display.rownames=TRUE)
        invisible(newdata)
        }   
    if (type=="resid") {
        xx       = x$data$value
        yy       = predict(x)-xx
        print(xx)
        print(yy)
        if (norm == "standardized") { yy = (yy-mean(yy, na.rm=TRUE))/sd(yy, na.rm=TRUE) }
        if (norm == "weighted")    { yy = yy/xx }
        my = max(abs(yy), na.rm=TRUE)*1.05
        print(yy)
        plot(yy~xx, pch=symbol, ylim=c(-my,my), log="x", xlab=paste(analyte, "concentration", units),
            ylab=paste("Residuals (", norm, ")", sep=""), ...) 
        abline(h=0, col="#88888899", lwd=1)
        }
}

predict.sigfit <- function(object, newdata=NULL, e.fit=FALSE, ...) {
    # Set stuff
    fit = object$fit
    analyte= object$analyte[1]
    model= object$model[1]
    type= object$model[2]
    
    # Determine what kind of data to use
    if (is.null(newdata)) { data    = object$data }
    else {
        if (any(class(newdata)=="ima") & ncol(newdata)>2) {
            if (object$file!=attr(newdata, "file")) stop("Can't predict data from a different run than the model was built on.")
            fulldata = newdata
            data  = as.data.frame(newdata[,paste(c("MFI","conc"), rep(analyte,2), sep=".")])
            names(data)[1:2] = c("MFI","value")
            }
        else {
            nn   = names(newdata)
            if (length(grep("MFI.",nn))>0) { 
                if (substr(nn[1],5,50)!=analyte) warning("Analyte name in model object does not match the data.") }
            data = newdata 
            }
        names(data)[1:2] = c("MFI","value")
        }

    # Calculations
    s       = coefficients(fit)
    if (toupper(substr(model,0,1))=="L") {
        if (type=="4")  z.preds = s[3] * ((s[2]/(data$MFI-s[1]))-1)^(1/s[4])
        if (type=="5")  z.preds = s[3] * ((s[2]/(data$MFI-s[1]))^(1/s[5])-1)^(1/s[4])
        if (e.fit) z.preds = data.frame(preds=z.preds, error = (z.preds - data$value) / data$value * 100)
        }
    if (toupper(substr(model,0,1))=="H") {
        if (type=="4")  suppressWarnings( { 
            z.preds = exp(s[3] - ((log10((s[2]-data$MFI)/(data$MFI-s[1])))/s[4])) } )
        if (type=="5")  suppressWarnings( { 
            z.preds = exp(s[3] - ((log10((((s[2]-s[1])/(data$MFI-s[1]))^(1/s[5]))-1))/s[4])) } )
        if (e.fit) z.preds = data.frame(preds=z.preds, error = (z.preds - data$value) / data$value * 100)
        }
        
    # Report back
    if (exists("fulldata")) {
        if (class(z.preds)=="data.frame") fulldata[,paste("pred", analyte, sep=".")] = z.preds[,1]
        else fulldata[,paste("pred", analyte, sep=".")] = z.preds
        if (e.fit) {
            IDx = unique(paste(fulldata$ID, fulldata$SPL))
            for (i1 in IDx) {
                bz = fulldata[paste(fulldata$ID, fulldata$SPL) %in% i1, paste("pred", analyte, sep=".")]
                fulldata[paste(fulldata$ID, fulldata$SPL) == i1, paste("pred", analyte, "cv", sep=".")] =
                    sd(bz, na.rm=TRUE)/mean(bz, na.rm=TRUE) * 100
                }  
            }
        invisible(fulldata)
        }
    else { return(z.preds) }
}

print.sigfit <- function(x, ...) {
    cat(ifelse(x$model["equation"]=="L","Logistic","Hill"),"type", 
    x$model["type"], "fit for", x$analyte["analyte"])
    if (!is.na(x$model[3])) cat(" with", x$model[3], "weighting.\n")
    else cat(".\n")
    cat("Data file:",x$file,"\n\n")
    print(x$fit)
    cat("\nBased on calibrators:\n")
    print(x$data)
    if (!all(is.na(x$stats[1]))) {
        cat("\nFit statistics for selected fits:\n")
        print(x$stats[,c(1:2,7:10)])
        }
}

summary.sigfit <- function(object, ...) {
    s = try(summary(object$fit), silent=TRUE)
    if (class(s)=="try-error") stop(ifelse(object$model["equation"]=="L","Logistic","Hill")," type ", 
    object$model["type"], " fit for ", object$model["analyte"]," - invalid.") 
    else {        
        d            = object$qcs
        if (nrow(d)>0) { d$weights    = d$use = NA }
        d            = rbind(object$data, d)
        d$predicted  = predict(object, newdata=d)
        dd           = aggregate(predicted~unlist(strsplit(rownames(d), split=" "))[seq(1,nrow(d)*2,2)], data=d, FUN=mesdci)
        dd           = dd[match(unique(unlist(strsplit(rownames(d), split=" "))[seq(1,nrow(d)*2,2)]),dd[,1]),]
        d$accuracy   = paste(formatC(((d$predicted-d$value)/d$value)*100, format="f", digits=2), "%", sep="")
        d$predicted  = formatC(d$predicted, format="f", digits=2)
        prec = NULL
        for (i1 in 1:nrow(dd)) { prec  = c(prec, paste(round(as.numeric(dd$predicted[i1,3]),2), "%", sep=""),"") }
        d$precision  = prec
        d.sts        = d[substr(rownames(d),0,2) == "St",]
        d.qcs        = d[substr(rownames(d),0,2) != "St",]
        
        cat(ifelse(object$model["equation"]=="L","Logistic","Hill"),"type", 
        object$model["type"], "fit for", object$analyte["analyte"])
        if (!is.na(object$model[3])) cat(" with", object$model[3], "weighting.")
        else cat(".")
        cat("\nData file:",object$file,"\n")
        print(s)
        cat("Calibrators accuracy & precision:\n")
        print(d.sts)
        if (nrow(d.qcs)>0) {
            cat("\nQCs accuracy & precision:\n")
            print(d.qcs)
            }
        if (!all(is.na(object$stats[1]))) {
            cat("\nFit statistics for selected fits:\n")
            print(object$stats)
            }
        invisible(list(fit=s, res=d, stats=object$stats))
        }
}

