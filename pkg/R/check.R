check <- function(x, ...) {
    z = class(x)[1]
    if (z %in% c("sigfit","ima")) {
        if (z=="sigfit") check.sigfit(x, ...)
        if (z=="ima") check.ima(x, ...)
    }
    else stop("Method \"check\" not defined for object of class: \"",z,"\".")
}

check.sigfit <- function(x) {
    analyte     = x$analyte["analyte"]
    d           = x$qcs
    if (nrow(d)>0) { d$weights   = d$use = NA }
    d           = rbind(x$data, d)
    d$Type      = "Standard"
    d$Type[substr(rownames(d),0,2)!="St"] = "QC"

    n.fit       = predict(x, newdata=d, e.fit=TRUE)
    n.fit$error = abs(n.fit$error)
    n.fit       = cbind(d, n.fit)
    n.fit$SPL   = unlist(strsplit(rownames(d), split=" "))[seq(1,nrow(d)*2,2)]
    cals        = unique(n.fit$SPL[n.fit$Type=="Standard" & n.fit$use==1])
    st.err      = c(median = median(abs(n.fit$error[n.fit$Type=="Standard"]), na.rm=TRUE),
    mean        = mean(abs(n.fit$error[n.fit$Type=="Standard"]), na.rm=TRUE))
    qc.err      = c(median = median(abs(n.fit$error[n.fit$Type=="QC"]), na.rm=TRUE),
                    mean   = mean(abs(n.fit$error[n.fit$Type=="QC"]), na.rm=TRUE))
    sse         = sum(((n.fit$error[n.fit$SPL %in% cals])*n.fit$use[n.fit$SPL %in% cals])^2, na.rm=TRUE)
    sse_2       = sum((n.fit$error[n.fit$SPL %in% cals[2:length(cals)-1]])^2, na.rm=TRUE)
    sse_3       = sum((n.fit$error[n.fit$SPL %in% cals[2:length(cals)-2]])^2, na.rm=TRUE)
    rsq         = as.vector(1 - sse / sum((n.fit$preds[n.fit$SPL %in% cals] - 
                  mean(n.fit$preds[n.fit$SPL %in% cals], na.rm=TRUE))^2, na.rm=TRUE))
    rsq_2       = as.vector(1 - sse_2 / sum((n.fit$preds[n.fit$SPL %in% cals[2:length(cals)-1]] - 
                  mean(n.fit$preds[n.fit$SPL %in% cals[2:length(cals)-1]], na.rm=TRUE))^2, na.rm=TRUE))
    rsq_3       = as.vector(1 - sse_3 / sum((n.fit$preds[n.fit$SPL %in% cals[2:length(cals)-2]] - 
                  mean(n.fit$preds[n.fit$SPL %in% cals[2:length(cals)-2]], na.rm=TRUE))^2, na.rm=TRUE))
    syx         = sqrt(deviance(x$fit)/(length(cals)-length(coef(x$fit))))
    rse         = sqrt(1/(nrow(n.fit[n.fit$Type=="Standard",])-length(coef(x$fit))) * sse)
    return(list(St.error=st.err, QC.error=qc.err, SSE=sse, sigma=rse, Syx=syx, r.squared=rsq))
}

check.ima <- function(x, analyte) {
    d           = as.data.frame(x[x$Type %in% c("Standard","QC"),])
    d           = d[, c("Type", "SPL", paste(c("RES","MFI","conc"), analyte, sep="."))]
    names(d)    = c("Type","SPL","preds","MFI","value")

    d$error     = abs((d$preds - d$value) / d$value * 100)

    n.fit       = d
    cals        = unique(n.fit$SPL[n.fit$Type=="Standard"])
    st.err      = c(median = median(abs(n.fit$error[n.fit$Type=="Standard"]), na.rm=TRUE),
                    mean   = mean(abs(n.fit$error[n.fit$Type=="Standard"]), na.rm=TRUE))
    qc.err      = c(median = median(abs(n.fit$error[n.fit$Type=="QC"]), na.rm=TRUE),
                    mean   = mean(abs(n.fit$error[n.fit$Type=="QC"]), na.rm=TRUE))
    sse         = sum(((n.fit$error[n.fit$SPL %in% cals]))^2, na.rm=TRUE)
    sse_2       = sum((n.fit$error[n.fit$SPL %in% cals[2:length(cals)-1]])^2, na.rm=TRUE)
    sse_3       = sum((n.fit$error[n.fit$SPL %in% cals[2:length(cals)-2]])^2, na.rm=TRUE)
    rsq         = as.vector(1 - sse / sum((n.fit$preds[n.fit$SPL %in% cals] - 
                  mean(n.fit$preds[n.fit$SPL %in% cals], na.rm=TRUE))^2, na.rm=TRUE))
    rsq_2       = as.vector(1 - sse_2 / sum((n.fit$preds[n.fit$SPL %in% cals[2:length(cals)-1]] - 
                  mean(n.fit$preds[n.fit$SPL %in% cals[2:length(cals)-1]], na.rm=TRUE))^2, na.rm=TRUE))
    rsq_3       = as.vector(1 - sse_3 / sum((n.fit$preds[n.fit$SPL %in% cals[2:length(cals)-2]] - 
                  mean(n.fit$preds[n.fit$SPL %in% cals[2:length(cals)-2]], na.rm=TRUE))^2, na.rm=TRUE))
    syx         = "not available"
    rse         = sqrt(1/(nrow(n.fit[n.fit$Type=="Standard",])- nrow(attr(a, "coefs"))) * sse)
    return(list(St.error=st.err, QC.error=qc.err, SSE=sse, sigma=rse, Syx=syx, r.squared=rsq))
}


