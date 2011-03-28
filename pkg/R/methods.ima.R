plot.ima <- function(data, type="cts", analyte=1, ref=0.25, cts.scale=350, ...) {
    # Basic definitions
    z.analytes = attr(data, "Analytes")
    cls  = c("#ff225577","#00ff3377","#ffaa1f77","#0000ff77","#2222cc77","#44444477","#44444477")
    i1   = analyte
    z.n  = as.numeric(attr(data, "Kit")[2])
    ord  = as.numeric(unlist(strsplit(as.character(data$Loc), split="(", fixed=T))[seq(1,nrow(data)*2,2)])

    # Define plot and reference matrices
    if (tolower(type)=="cts") {
        mtx = matrix(data[order(ord), paste("CTS",z.analytes[i1], sep=".")], nrow=8, ncol=12, byrow=F)
        mtx = mtx/cts.scale
        mtl = paste("Counts of", z.analytes[i1])
        xl  = "Reference circle: 20 counts"
        }   
    if (tolower(type)=="cv") {
        if (!is.na(match(paste("pred", z.analytes[i1], sep="."), names(data)))) {
            if (is.na(match(paste("pred", z.analytes[i1], "cv", sep="."), names(data)))) {
                IDx = unique(paste(data$ID, data$SPL))
                for (i2 in IDx) {
                    bz = data[paste(data$ID, data$SPL) %in% i2, paste("pred", z.analytes[i1], sep=".")]
                    data[paste(data$ID, data$SPL) == i2, paste("pred", z.analytes[i1], "cv", sep=".")] =
                        sd(bz, na.rm=T)/mean(bz, na.rm=T) * 100
                    }  
                }
            mtx = matrix(data[order(ord), paste("pred", z.analytes[i1], "cv", sep=".")], ncol=12, byrow=F)
            mtl = paste("Precision of results for", z.analytes[i1])
            }
        else {
            CV  = aggregate(data[,paste("MFI",z.analytes[i1], sep=".")], 
                  by = list(ID=data$ID,Type=data$Type,SPL=data$SPL), FUN=cv)
            CV  = merge(data, CV,  by=c("ID","Type","SPL"), all=T)
            mtx = matrix(CV[order(ord), "x"], ncol=12, byrow=F)
            mtl = paste("Precision of MFI for", z.analytes[i1])
            }
        mtx = mtx/100
        xl  = paste("Reference circle: ", ref*100,"% CV", sep="")
        }
    if (tolower(type)=="accuracy") {
        if (!is.na(match(paste("pred", z.analytes[i1], sep="."), names(data)))) {
            mtx = matrix(abs(1 - data[order(ord), paste("pred", z.analytes[i1], sep=".")] / 
            data[order(ord), paste("conc", z.analytes[i1], sep=".")]), ncol=12, byrow=F)
            mtl = paste("Accuracy of results for", z.analytes[i1])
            }
        else stop("Accuracy cannot be derived from an object without predictions")
        xl  = paste("Reference circle: ", ref*100,"% deviation from set value", sep="")
        }
    if (z.n<96) { mtx[(z.n+1):96] = NA }
    if (ref<=0) { xl="" }
    mtx[tolower(data[order(ord),"SPL"])=="background"] = NA
    suppressWarnings( { mtx = sqrt(mtx/pi) } )
    
    # Make plot
    {
    plot(0,0, type="n", xlab="", ylab="", xlim=c(1,12), ylim=c(0.75,8.25), axes=FALSE, ...)
    mtext(xl, side=1, line=0, adj=0.5)
    symbols(rep(c(1:12), each=8), rep(8:1, 12), circles=as.vector(mtx), inches=F, 
        bg=cls[as.numeric(data[order(ord),"Type"])], ylab="", xlab="", add=T, xpd=NA)
    axis(3, at=1:12, labels=1:12, line=-0.5, lwd=0, lwd.ticks=0, col.ticks="gray")
    axis(2, at=1:8, labels=LETTERS[8:1], line=-0.5, lwd=0, lwd.ticks=0, col.ticks="gray", las=1)
    abline(h=c(1:7)+0.5, v=c(1:11)+0.5, col="#55555555", lty=2)
    title(main=mtl, line=2.5)
    }
    
    # Plot reference circles
    if (ref>0) {
        if (tolower(type)=="cts") { mtx = sqrt(matrix(rep(20/cts.scale,96), ncol=12, byrow=F)/pi) }
        else { mtx = sqrt(matrix(rep(ref,96), ncol=12, byrow=F)/pi) }
        if (z.n<96) { mtx[(z.n+1):96] = NA }
        symbols(rep(c(1:12), each=8), rep(8:1, 12), circles=as.vector(mtx), inches=F, 
            fg="#55555555", ylab="", xlab="", add=T, xpd=NA)
        }
    if (!is.na(match(paste("pred", z.analytes[i1], sep="."), names(data)))) {
        if (tolower(type)=="cv") { 
            mtx = matrix(data[order(ord), paste("pred", z.analytes[i1], "cv", sep=".")], ncol=12, byrow=F) 
            mtx[!is.na(mtx) & mtx>ref*100] = NA
            mtext("\"X\" - exceeding the limit or missing" , side=1, line=1, adj=0.5, cex=0.65)
            }
        if (tolower(type)=="accuracy") { 
            mtx = matrix(data[order(ord), paste("conc", z.analytes[i1], sep=".")], ncol=12, byrow=F) 
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

print.ima <- function(x) {
    cat("Immunoassay run file:", attr(x, "file"), "\n")
    cat("Number of samples:", as.numeric(attr(x, "Kit")[2]), "\n")
    cat("Kit Lot Number:", attr(x, "Kit")[1], "\n")
    cat("Analytes:", attr(x,"Analytes"), "\n")
    cat("Background:", attr(x,"Background"), "\n")
    cat("Date processed:", attr(x,"Date"), "\n")
    cat("On:", attr(x,"Assay")[1], "  S/N:", attr(x,"Assay")[2], "\n")
    cat("By:", attr(x,"Operator"), "\n\n")
    d = as.data.frame(x)
    rownames(d) = paste(d$SPL, " (", unlist(strsplit(as.character(d$Loc), split="(", fixed=T))[seq(2,nrow(d)*2,2)], sep="")
    print(d[,-match(c("SPL","Loc","ID"), names(d))])
}


