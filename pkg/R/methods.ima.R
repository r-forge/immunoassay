plot.ima    <- function(x, what="mfi", type="precision", analyte=1, ref=0.25, cts.scale=350, ...) {
    if (!(tolower(what) %in% c("cts","mfi","res","pred"))) stop("Unsupported result type.")
    if (!(tolower(type) %in% c("precision","accuracy")))  stop("Unsupported plot type.")
 
    # Utility functions
    CV <- function(dat, column, ids) {
        CV  = aggregate(dat[,column], by = ids, FUN= function(x) { sd(x, na.rm=T)/mean(x, na.rm=T)*100 } )
        return(CV)
    }
    
    # Basic definitions
    z.analytes = attr(x, "Analytes")
    cls  = c("#ff225577","#00ff3377","#ffaa1f77","#0000ff77","#2222cc77","#44444477","#44444477")
    i1   = analyte
    z.n  = as.numeric(attr(x, "Lot")[2])
    ord  = as.numeric(matrix(unlist(strsplit(as.character(x$Loc),"(", fixed=T)), ncol=2, byrow=T)[,1])
    
    # Error checks
    if (z.n>96) stop("Unsupported plate layout")
    if (what=="res" & length(grep(paste("RES.", z.analytes[i1], sep=""), names(x)))!=1) stop("Results are missing for analyte ", z.analytes[i1])
    if (what=="pred" & length(grep(paste("pred.", z.analytes[i1], sep=""), names(x)))!=1) stop("Predictions are missing for analyte ", z.analytes[i1])
    if (what=="mfi" & type!="precision") stop("MFI plot can only be of type \"precision\"")
    
    # Conversions
    if (tolower(what)!="pred") what=toupper(what)
    if (what=="CTS") type=""
    
    # Define plot and reference matrices
    if (what=="CTS") {
        mtx = matrix(x[order(ord), paste(what, z.analytes[i1], sep=".")], nrow=8, ncol=12, byrow=FALSE)
        mtx = mtx/cts.scale
        mtl = paste("Counts of", z.analytes[i1])
        xl  = "Reference circle: 50 counts"
    } else {
        if (tolower(type)=="precision") {
            if (!is.na(match(paste(what, z.analytes[i1], sep="."), names(x)))) {
                if (is.na(match(paste(what, z.analytes[i1], "cv", sep="."), names(x)))) {
                    bz = CV(x, paste(what, z.analytes[i1], sep="."), list(x$SPL))
                    names(bz)[2] = paste(what, z.analytes[i1], "cv", sep=".")
                    x  = merge(x, bz, by.x="SPL", by.y="Group.1", all.x=T)
                }
                ord  = as.numeric(matrix(unlist(strsplit(as.character(x$Loc),"(", fixed=T)), ncol=2, byrow=T)[,1])
                mtx = matrix(x[order(ord), paste(what, z.analytes[i1], "cv", sep=".")], ncol=12, byrow=FALSE)
                mtl = paste("Precision of", ifelse(what=="RES","results", ifelse(what=="pred", "predictions", what) ),"for", z.analytes[i1])
            }
            mtx = mtx/100
            xl  = paste("Reference circle: ", ref*100,"% CV", sep="")
        }
        if (tolower(type)=="accuracy") {
            if (!is.na(match(paste(what, z.analytes[i1], sep="."), names(x)))) {
                mtx = matrix(abs(1 - x[order(ord), paste(what, z.analytes[i1], sep=".")] / 
                x[order(ord), paste("conc", z.analytes[i1], sep=".")]), ncol=12, byrow=FALSE)
                mtl = paste("Accuracy of", ifelse(what=="RES","results", "predictions"),"for", z.analytes[i1])
                }
            xl  = paste("Reference circle: ", ref*100,"% deviation from nominal value", sep="")
        }
    }
    if (z.n<96) { mtx[(z.n+1):96] = NA }
    if (ref<=0) { xl="" }
    mtx[tolower(x[order(ord),"SPL"])=="background"] = NA
    suppressWarnings( { mtx = sqrt(mtx/pi) } )
    
    # Make plot
    {
    plot(0,0, type="n", xlab="", ylab="", xlim=c(1,12), ylim=c(0.75,8.25), axes=FALSE, ...)
    mtext(xl, side=1, line=0.3, adj=0.5)
    symbols(rep(c(1:12), each=8), rep(8:1, 12), circles=as.vector(mtx), inches=FALSE, 
        bg=cls[as.numeric(x[order(ord),"Type"])], ylab="", xlab="", add=TRUE, xpd=NA)
    axis(3, at=1:12, labels=1:12, line=-0.5, lwd=0, lwd.ticks=0, col.ticks="gray")
    axis(2, at=1:8, labels=LETTERS[8:1], line=-0.5, lwd=0, lwd.ticks=0, col.ticks="gray", las=1)
    abline(h=c(1:7)+0.5, v=c(1:11)+0.5, col="#55555555", lty=2)
    title(main=mtl, line=2.5)
    }
    
    # Plot reference circles
    if (ref>0) {
        if (what=="CTS") { mtx = sqrt(matrix(rep(50/cts.scale,96), ncol=12, byrow=FALSE)/pi) }
        else { mtx = sqrt(matrix(rep(ref,96), ncol=12, byrow=FALSE)/pi) }
        if (z.n<96) { mtx[(z.n+1):96] = NA }
        symbols(rep(c(1:12), each=8), rep(8:1, 12), circles=as.vector(mtx), inches=FALSE, 
            fg="#55555555", ylab="", xlab="", add=TRUE, xpd=NA)
    }
    if (what!="cts" & !is.na(match(paste(what, z.analytes[i1], sep="."), names(x)))) {
        if (tolower(type)=="precision") { 
            mtx = matrix(x[order(ord), paste(what, z.analytes[i1], "cv", sep=".")], ncol=12, byrow=FALSE) 
            mtx[!is.na(mtx) & mtx>ref*100] = NA
            mtext("\"X\" - exceeding the limit or missing" , side=1, line=1, adj=0.5, cex=0.65)
        }
        if (tolower(type)=="accuracy") { 
            mtx = matrix(x[order(ord), paste("conc", z.analytes[i1], sep=".")], ncol=12, byrow=FALSE) 
            mtx[c(1,9)] = NA
            mtext("\"X\" - not evaluated" , side=1, line=1, adj=0.5, cex=0.65)
        }
        mtx[!is.na(mtx)] = 0
        mtx[is.na(mtx)]  = 4
        mtx[mtx==0]  = NA 
        if (z.n<96) { mtx[(z.n+1):96] = NA }
        points(rep(c(1:12), each=8), rep(8:1, 12), pch=as.vector(mtx), cex=2)
    }
}

print.ima   <- function(x, ...) {
    cat("Immunoassay run file:", attr(x, "file"), "\n")
    cat("Immunoassay Session/Batch:", attr(x, "Assay")[3], "\n")
    cat("Number of samples:", as.numeric(attr(x, "Lot")[2]), "\n")
    cat("Reagents Lot Number:", attr(x, "Lot")[1], "\n")
    cat("Analytes:", attr(x,"Analytes"), "\n")
    cat("Background:", attr(x,"Background"), "\n")
    cat("Date processed:", attr(x,"Date"), "\n")
    cat("On:", attr(x,"Assay")[1], "  S/N:", attr(x,"Assay")[2], "\n")
    cat("By:", attr(x,"Operator"), "\n\n")
    d = as.data.frame(x)
    rownames(d) = paste(d$SPL, " (", unlist(strsplit(as.character(d$Loc), split="(", fixed=TRUE))[seq(2,nrow(d)*2,2)], sep="")
    print(d[,-match(c("SPL","Loc","ID"), names(d))], ...)
}

predict.ima <- function(object, analyte=1, newdata=NULL, ...) {
    if (length(grep("xPONENT", attr(object, "Assay")[1], fixed=T))!=1) stop("Fit data not available")
    if (analyte>length(attr(object,"Analytes"))) stop("Analyte does not exist")
    
    inv <- function(y, a, b, c, d, f) {
        c * (((b-a)/(y-a))^(1/f)-1)^(1/d)
    }

    an      = attr(object, "Analytes")[analyte]
    coefs   = as.numeric(t(attr(object, "coefs"))[an,])
    
    if (is.null(newdata)) {
        preds = inv(object[,paste("MFI",an, sep=".")], a=coefs[1], b=coefs[2], c=coefs[3], d=coefs[4], f=coefs[5])
    } else {
        preds = inv(newdata, a=coefs[1], b=coefs[2], c=coefs[3], d=coefs[4], f=coefs[5])
    }
    
    return(preds)
}

summary.ima <- function(object, analyte="all", result="res", type="fit", ...) {
    if (!(tolower(result) %in% c("res","pred"))) stop("Unsupported result type")
    if (!(tolower(type)   %in% c("fit","data"))) stop("Unsupported summary type")
    if (!is.numeric(analyte) & analyte!="all") stop("Incorrect analyte selected.")
    
    cat("Immunoassay run file:", attr(object, "file"), "\n")
    cat("Immunoassay Session/Batch:", attr(object, "Assay")[3], "\n")
    cat("Number of samples:", as.numeric(attr(object, "Lot")[2]), "\n")
    cat("Reagents Lot Number:", attr(object, "Lot")[1], "\n")
    cat("Analytes:", attr(object,"Analytes"), "\n")
    cat("Background:", attr(object,"Background"), "\n")
    cat("Date processed:", attr(object,"Date"), "\n")
    cat("On:", attr(object,"Assay")[1], "  S/N:", attr(object,"Assay")[2], "\n")
    cat("By:", attr(object,"Operator"), "\n\n")

    if (!is.numeric(analyte) & analyte=="all") an = attr(object, "Analytes")
    else an = attr(object, "Analytes")[analyte]
    
    if (type=="fit") {
        if (length(grep("xPONENT", attr(object, "Assay")[1], fixed=T))==1) { 
            a1 = unlist(attr(object, "fit")$fit)
            
            tab1 = cbind(levels(a1)[a1], unclass(attr(object, "fit")$r2))
            rownames(tab1) = an
            colnames(tab1) = c("Fit type","R squared")
            cat("Fitting summary:\n")
            print(tab1, quote=FALSE)
            tab2 = t(attr(object, "coefs"))[(1:length(an)+1),]
            colnames(tab2) = c("a","b","c","d","f")
            cat("\nCoefficients:\n")
            print(tab2, quote="FALSE")
            cat("\n")
        } else {
            cat("\nFit data not available\n\n")
        }

        if (tolower(result)=="res") result="RES" else result="pred"
        d = as.data.frame(object[object$Type %in% c("Standard","QC"), c("SPL","Loc",paste(rep(c("MFI", result, "conc"), each=length(an)), an, sep="."))])
        d$Loc = paste("(",matrix(unlist(strsplit(as.character(d$Loc), ",", fixed=T)), ncol=2, byrow=T)[,2], sep="")

        for (i1 in an) {          
            dd   = d[, paste(c("MFI","conc",result), i1, sep=".")]
            colnames(dd) = c("MFI","value","result")
            rownames(dd) = paste(d$SPL,d$Loc)
            dd$accuracy  = paste(formatC((dd$result-dd$value)/dd$value * 100, format="f", digits=2), "%", sep="")
            dd$result    = round(dd$result, 2)
            pr           = aggregate(dd$result, list(d$SPL), FUN = function(x) { sd(x, na.rm=T)/mean(x, na.rm=T)*100 })
            dd$precision = ""
            dd[match(pr$Group.1, d$SPL), "precision"] = paste(formatC(pr$x, format="f", digits=2), "%", sep="")

            cat("Accuracy & precision of Standards and QCs for ", i1,":\n", sep="")
            print(dd)
            cat("\n")
        }
    } else {
        summary(as.data.frame(object))
    }
}
