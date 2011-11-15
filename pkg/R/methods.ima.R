plot.ima <- function(x, type="cts", analyte=1, ref=0.25, cts.scale=350, ...) {
    # Basic definitions
    z.analytes = attr(x, "Analytes")
    cls  = c("#ff225577","#00ff3377","#ffaa1f77","#0000ff77","#2222cc77","#44444477","#44444477")
    i1   = analyte
    z.n  = as.numeric(attr(x, "Lot")[2])
    ord  = as.numeric(unlist(strsplit(as.character(x$Loc), split="(", fixed=TRUE))[seq(1,nrow(x)*2,2)])

    # Define plot and reference matrices
    if (tolower(type)=="cts") {
        mtx = matrix(x[order(ord), paste("CTS",z.analytes[i1], sep=".")], nrow=8, ncol=12, byrow=FALSE)
        mtx = mtx/cts.scale
        mtl = paste("Counts of", z.analytes[i1])
        xl  = "Reference circle: 20 counts"
        }   
    if (tolower(type)=="cv") {
        if (!is.na(match(paste("pred", z.analytes[i1], sep="."), names(x)))) {
            if (is.na(match(paste("pred", z.analytes[i1], "cv", sep="."), names(x)))) {
                IDx = unique(paste(x$ID, x$SPL))
                for (i2 in IDx) {
                    bz = x[paste(x$ID, x$SPL) %in% i2, paste("pred", z.analytes[i1], sep=".")]
                    x[paste(x$ID, x$SPL) == i2, paste("pred", z.analytes[i1], "cv", sep=".")] =
                        sd(bz, na.rm=TRUE)/mean(bz, na.rm=TRUE) * 100
                    }  
                }
            mtx = matrix(x[order(ord), paste("pred", z.analytes[i1], "cv", sep=".")], ncol=12, byrow=FALSE)
            mtl = paste("Precision of results for", z.analytes[i1])
            }
        else {
            CV  = aggregate(x[,paste("MFI",z.analytes[i1], sep=".")], 
                  by = list(ID=x$ID,Type=x$Type,SPL=x$SPL), FUN=cv)
            CV  = merge(x, CV,  by=c("ID","Type","SPL"), all=TRUE)
            mtx = matrix(CV[order(ord), "x"], ncol=12, byrow=FALSE)
            mtl = paste("Precision of MFI for", z.analytes[i1])
            }
        mtx = mtx/100
        xl  = paste("Reference circle: ", ref*100,"% CV", sep="")
        }
    if (tolower(type)=="accuracy") {
        if (!is.na(match(paste("pred", z.analytes[i1], sep="."), names(x)))) {
            mtx = matrix(abs(1 - x[order(ord), paste("pred", z.analytes[i1], sep=".")] / 
            x[order(ord), paste("conc", z.analytes[i1], sep=".")]), ncol=12, byrow=FALSE)
            mtl = paste("Accuracy of results for", z.analytes[i1])
            }
        else stop("Accuracy cannot be derived from an object without predictions")
        xl  = paste("Reference circle: ", ref*100,"% deviation from set value", sep="")
        }
    if (z.n<96) { mtx[(z.n+1):96] = NA }
    if (ref<=0) { xl="" }
    mtx[tolower(x[order(ord),"SPL"])=="background"] = NA
    suppressWarnings( { mtx = sqrt(mtx/pi) } )
    
    # Make plot
    {
    plot(0,0, type="n", xlab="", ylab="", xlim=c(1,12), ylim=c(0.75,8.25), axes=FALSE, ...)
    mtext(xl, side=1, line=0, adj=0.5)
    symbols(rep(c(1:12), each=8), rep(8:1, 12), circles=as.vector(mtx), inches=FALSE, 
        bg=cls[as.numeric(x[order(ord),"Type"])], ylab="", xlab="", add=TRUE, xpd=NA)
    axis(3, at=1:12, labels=1:12, line=-0.5, lwd=0, lwd.ticks=0, col.ticks="gray")
    axis(2, at=1:8, labels=LETTERS[8:1], line=-0.5, lwd=0, lwd.ticks=0, col.ticks="gray", las=1)
    abline(h=c(1:7)+0.5, v=c(1:11)+0.5, col="#55555555", lty=2)
    title(main=mtl, line=2.5)
    }
    
    # Plot reference circles
    if (ref>0) {
        if (tolower(type)=="cts") { mtx = sqrt(matrix(rep(20/cts.scale,96), ncol=12, byrow=FALSE)/pi) }
        else { mtx = sqrt(matrix(rep(ref,96), ncol=12, byrow=FALSE)/pi) }
        if (z.n<96) { mtx[(z.n+1):96] = NA }
        symbols(rep(c(1:12), each=8), rep(8:1, 12), circles=as.vector(mtx), inches=FALSE, 
            fg="#55555555", ylab="", xlab="", add=TRUE, xpd=NA)
        }
    if (!is.na(match(paste("pred", z.analytes[i1], sep="."), names(x)))) {
        if (tolower(type)=="cv") { 
            mtx = matrix(x[order(ord), paste("pred", z.analytes[i1], "cv", sep=".")], ncol=12, byrow=FALSE) 
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

print.ima <- function(x, ...) {
    cat("Immunoassay run file:", attr(x, "file"), "\n")
    cat("Number of samples:", as.numeric(attr(x, "Lot")[2]), "\n")
    cat("Reagents Lot Number:", attr(x, "Lot")[1], "\n")
    cat("Analytes:", attr(x,"Analytes"), "\n")
    cat("Background:", attr(x,"Background"), "\n")
    cat("Date processed:", attr(x,"Date"), "\n")
    cat("On:", attr(x,"Assay")[1], "  S/N:", attr(x,"Assay")[2], "\n")
    cat("By:", attr(x,"Operator"), "\n\n")
    d = as.data.frame(x)
    rownames(d) = paste(d$SPL, " (", unlist(strsplit(as.character(d$Loc), split="(", fixed=TRUE))[seq(2,nrow(d)*2,2)], sep="")
    print(d[,-match(c("SPL","Loc","ID"), names(d))])
}


