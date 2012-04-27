read.kits <- function(path, file) {
    if (!file.exists(path=paste(path, file, sep="/"))) stop(paste("No \"",file,"\" file found in the \"",path,"\" path.", sep=""))
    immunoassay.kits      = read.csv(paste(path, file, sep="/"))
    immunoassay.kits$ID   = factor(paste(immunoassay.kits$Lot, immunoassay.kits$Analyte))
    for (i1 in levels(immunoassay.kits$ID)) {
        immunoassay.kits$ConA[immunoassay.kits$ID==i1] = mean(unlist(immunoassay.kits[immunoassay.kits$ID==i1, c("ConA.lo","ConA.hi")]))
        immunoassay.kits$ConB[immunoassay.kits$ID==i1] = mean(unlist(immunoassay.kits[immunoassay.kits$ID==i1, c("ConB.lo","ConB.hi")]))
        }
    immunoassay.kits$ID = NULL
    return(immunoassay.kits)
}


read.multiplex <- function(path, file, analytes="all") {

    # Function definitions    
    extract.xP <- function(filename) {
        # Read run info
        a= scan(filename, what="character", sep=",", quiet=TRUE)
        z.n= as.numeric(a[grep("Samples",a)[1]+1])*2
        cnt         = 1
        repeat {
            z.n     = as.numeric(a[grep("Samples",a)[cnt]+1])*2
            if (!is.na(z.n)) break;
            cnt     = cnt+1
        }

        z.sn= a[grep("SN",a)[1]+1]
        z.platform  = paste(a[grep("xPONENT",a)[1]], a[grep("xPONENT",a)[1]+2])
        z.date= a[grep("BatchStopTime",a)[1]+1]
        z.operator  = a[grep("Operator",a)[1]+1]
        z.session   = a[grep("Batch",a)[1]+1]
        z.lot       = as.character(a[grep("Standard1", a)[1]+2])
        
        # Find locations of RESULTS, MFI, COUNTS and Calibrators tables
        for (i1 in 1:5000){
            a = scan(filename, what="character", sep=",", skip=i1, nlines=1, quiet=TRUE)
            if (length(a)>0) {
                if (a[1]=="DataType:") {
                    if (a[2]=="Median") { z.MFI.skip=i1+1; next }
                    if (a[2]=="Count")  { z.CTS.skip=i1+1; next }
                    if (a[2]=="Result") { z.RES.skip=i1+1; next }
                    if (a[2]=="Units")  { z.UNI.skip=i1+1; next }
                    if (a[2]=="Standard Expected Concentration") { z.CAL.skip=i1+1; next }
                    if (a[2]=="Control Range - Low")   { z.QCL.skip=i1+1; next }
                    if (a[2]=="Control Range - High")  { z.QCH.skip=i1+1; next }
                    if (a[2]=="Analysis Coefficients") { z.COE.skip=i1+1; next }
                    if (a[2]=="Analysis Types")        { z.FIT.skip=i1+1; next }
                    if (a[2]=="R^2")                   { z.R2.skip=i1+1; break }
                }
            }
        }
            
        # Find analyte names
        ll = scan(filename, skip=z.MFI.skip, nlines=1, what="character", sep=",", quiet=TRUE)
        analytes = ll[(match("Sample",ll)+1):(match("Total Events",ll)-1)]
        
        # Finish
        return(list(n=z.n,lot=z.lot,sn=z.sn,ses=z.session,platform=z.platform,date=z.date,operator=z.operator,
            skip=c(MFI=z.MFI.skip, CTS=z.CTS.skip, RES=z.RES.skip, UNI=z.UNI.skip, CAL=z.CAL.skip, 
            QCL=z.QCL.skip, QCH=z.QCH.skip, COEF=z.COE.skip, FIT=z.FIT.skip, R2=z.R2.skip), analytes=analytes))
    }

    extract.is <- function(filename) {
        # Read run info
        a= scan(filename, what="character", sep=",", quiet=TRUE)
        
        # Word "Samples" count
        cnt         = 1
        repeat {
            z.n     = as.numeric(a[grep("Samples",a)[cnt]+1])
            if (!is.na(z.n)) break;
            cnt     = cnt+1
        }
        z.sn        = a[grep("SN",a)[1]+1]
        z.platform  = a[grep("Program",a)[1]+1]
        z.date= a[grep("BatchStopTime",a)[1]+1]
        z.operator  = a[grep("Operator",a)[1]+1]
        z.session   = a[grep("Session",a)[1]+1]
        z.lot       = as.character(a[grep("Assay Standard", a)[1]+1])
        
        # Find locations of MFI and COUNTS tables
        for (i1 in 1:1000){
            a = scan(filename, what="character", sep=",", skip=i1, nlines=1, quiet=TRUE)
            if (length(a)>0) {
                if (a[1]=="DataType:" & a[2]=="Median") { z.MFI.skip=i1+1 }
                if (a[1]=="DataType:" & a[2]=="Result") { z.RES.skip=i1+1 }
                if (a[1]=="DataType:" & a[2]=="Count")  { z.CTS.skip=i1+1; break }
                }
            }

        # Find analyte names
        ll = scan(filename, skip=z.MFI.skip, nlines=1, what="character", sep=",", quiet=TRUE)
        analytes = ll[(match("Sample",ll)+1):(match("Total Events",ll)-1)]
        
        # Finish
        return(list(n=z.n,lot=z.lot,sn=z.sn,ses=z.session,platform=z.platform,date=z.date,operator=z.operator,
            skip=c(MFI=z.MFI.skip, CTS=z.CTS.skip, RES=z.RES.skip), analytes=analytes))
    }

    # Process path and file name
    if (tolower(substr(file, nchar(file)-3, nchar(file)))==".csv") {
        filename = paste(path, file, sep="")
        file     = substr(file, 0, nchar(file)-4)
    }
    else filename = paste(path, file, ".csv", sep="")
    if (!file.exists(filename)) stop(paste("No \"",file,"\" file found in the \"", path, "\" path.", sep=""))

    # Read first row of data
    a = scan(filename, what="character", sep=",", nlines=1, quiet=TRUE)
        
    # See if generated by IS or xPONENT software, and get the run info
    xP = length(grep("xPONENT",a)) 
    if (xP == 0 & length(grep("IS",a))==0) stop("Not a supported multiplex run file.")
    if (xP > 0) z.info = extract.xP(filename) 
    else z.info = extract.is(filename)
    
    # Set-up analyte names
    if (all(analytes=="all")) analytes = z.info$analytes
    else analytes = z.info$analytes[analytes]

    # Get analyte data & QC ranges
    if (xP==1) {
        # read calibrators and convert to common structure
        z.CAL = read.table(filename, sep=",", fill=TRUE, skip=z.info$skip["CAL"], nrows=12, header=TRUE)
        z.CAL = z.CAL[1:(match("DataType:", z.CAL[,1])-1),]
        cl.r  = data.frame(t(z.CAL)[1:length(z.info$analytes)+1,])
        names(cl.r) = paste("St", 1:ncol(cl.r), sep="")
        for (i1 in 1:ncol(cl.r)) cl.r[,i1] =  as.numeric(levels(cl.r[,i1])[cl.r[,i1]])
        cl.r$Analyte = z.info$analytes
        # read QCs and convert to common structure
        z.QCL = read.table(filename, sep=",", fill=TRUE, skip=z.info$skip["QCL"], nrows=2, header=TRUE)
        z.QCH = read.table(filename, sep=",", fill=TRUE, skip=z.info$skip["QCH"], nrows=2, header=TRUE)
        qc.r = data.frame(cbind(t(z.QCL[2:ncol(z.QCL)]), t(z.QCH[2:ncol(z.QCH)])))
        names(qc.r) = c("ConA.lo","ConB.lo","ConA.hi","ConB.hi"); 
        qc.r$Analyte = z.info$analytes
        for (i1 in qc.r$Analyte) {
            qc.r$ConA[qc.r$Analyte==i1] = mean(unlist(qc.r[qc.r$Analyte==i1, c("ConA.lo","ConA.hi")]))
            qc.r$ConB[qc.r$Analyte==i1] = mean(unlist(qc.r[qc.r$Analyte==i1, c("ConB.lo","ConB.hi")]))
        }
        # read coefficients
        z.COEF = read.table(filename, sep=",", fill=TRUE, skip=z.info$skip["COEF"], nrows=10, header=TRUE)
        z.COEF = z.COEF[1:(match("DataType:", z.COEF[,1])-1),]
        names(z.COEF) = c("Coefficient", z.info$analytes)
        # read units
        units = as.character(unlist(read.table(filename, sep=",", fill=TRUE, skip=z.info$skip["UNI"], 
            nrows=2, header=TRUE)[2,(1:length(analytes))+1]))
        cl.r = merge(cl.r, qc.r, by="Analyte")
        # read fit stats
        z.FIT = list(fit = read.table(filename, sep=",", fill=TRUE, skip=z.info$skip["FIT"], nrows=1, header=TRUE)[,2:(length(analytes)+1)],
                     r2  = read.table(filename, sep=",", fill=TRUE, skip=z.info$skip["R2"], nrows=1, header=TRUE)[,2:(length(analytes)+1)])
    } else {
        if (!exists("immunoassay.kits")) stop(paste("No \"immunoassay.kits\" data.frame found.", sep=""))
        kit  = immunoassay.kits[immunoassay.kits$Lot==as.character(z.info$lot),]
        if (nrow(kit)>0) { 
            if (!(all(analytes %in% kit$Analyte))) stop("Analytes from file \"", file, "\" not defined in \"immunoassay.kits\" data.frame.", sep="") 
            units = as.character(kit$Unit[as.character(kit$Lot)==z.info$lot & kit$Analyte %in% analytes]) 
            qc.r  = as.data.frame(kit[as.character(kit$Lot)==z.info$lot & kit$Analyte %in% analytes, 
                c("Analyte","ConA.lo","ConA.hi","ConB.lo","ConB.hi")])
            cl.r  = kit[as.character(kit$Lot)==z.info$lot & (kit$Analyte %in% analytes), ]
            z.COEF = NULL
            z.FIT  = NULL
            }
        else stop("Lot # not defined in \"immunoassay.kits\" data.frame.")
    }

    # Read in the tables
    z.MFI = read.table(filename, sep=",", fill=TRUE, skip=z.info$skip["MFI"], nrows=z.info$n, header=TRUE)
    z.MFI = z.MFI[,1:(ncol(z.MFI)-1)]
    z.CTS = read.table(filename, sep=",", fill=TRUE, skip=z.info$skip["CTS"], nrows=z.info$n, header=TRUE)
    z.CTS = z.CTS[,1:(ncol(z.CTS)-1)]
    z.RES = read.table(filename, sep=",", fill=TRUE, skip=z.info$skip["RES"], nrows=z.info$n, header=TRUE)
    z.RES = z.RES[,1:(ncol(z.RES)-1)]
    
    # Get rid of extra columns & merge
    names(z.MFI) = c("Loc","SPL",paste("MFI",analytes,sep="."))
    names(z.CTS) = c("Loc","SPL",paste("CTS",analytes,sep="."))
    names(z.RES) = c("Loc","SPL",paste("RES",analytes,sep="."))
    z.MFI        = z.MFI[,c("Loc","SPL",paste("MFI",analytes,sep="."))]
    z.CTS        = z.CTS[,c("Loc","SPL",paste("CTS",analytes,sep="."))]
    z.RES        = z.RES[,c("Loc","SPL",paste("RES",analytes,sep="."))]
    
    # Merge & Convert all to numeric
    z.run        = merge(z.MFI,z.CTS, by=c("Loc","SPL"))
    z.run        = merge(z.run,z.RES, by=c("Loc","SPL"))
    for (i1 in 1:(length(analytes)*3)) {
        if (class(z.run[,i1+2])!="numeric") {
            z.t          = unlist(strsplit(tolower(as.character(z.run[,i1+2])), " pg/ml", fixed=T))
            if (length(z.t)<z.info$n) z.t = c(z.t, rep(NA, z.info$n-length(z.t)))
            z.run[,i1+2] = as.numeric(z.t)
        }
    }

    # Add SPL class factor & re-create SPL names for retest samples
    z.run$SPL = as.character(z.run$SPL)
    z.run$SPL   = gsub(pattern = "S([A-Z])", replacement = "S\\L\\1", x=z.run$SPL, perl=TRUE)  # for "ST1" standard names
    z.run$Type  = "Sample"
    for (i1 in 1:nrow(z.run)) {
        if (substr(z.run$SPL[i1], nchar(z.run$SPL[i1]), nchar(z.run$SPL[i1])+1) == "R") {
            z.run$Type[i1] = "Retest"
            z.run$SPL[i1]  = substr(z.run$SPL[i1], 0, nchar(z.run$SPL[i1])-1)
            }
        }
    z.run$Type[substr(z.run$SPL,0,2)=="Co"] = "QC"
    z.run$Type[substr(z.run$SPL,0,1)=="Q"]  = "Pool"
    z.run$Type[tolower(substr(z.run$SPL,0,4))=="pool"]  = "Pool"
    z.run$Type[substr(z.run$SPL,0,2)=="St"] = "Standard"
    z.run$Type[substr(z.run$SPL,0,2) %in% c("Ba","ba")] = "Background"
    z.run$Type[z.run$SPL=="0" | tolower(substr(z.run$SPL,0,2))=="no"]     = "Missing"
    z.run$Type = ordered(z.run$Type, levels=c("Background","Standard","QC","Pool","Sample","Retest","Missing","Eliminated"))

    # Sort columns & rows
    z.run$ID = file
    z.run    = z.run[order(z.run$Type, z.run$SPL, as.numeric(unlist(strsplit(as.character(z.run$Loc), 
                  split="(", fixed=TRUE))[seq(1,nrow(z.run)*2,2)])), c("ID","Type","Loc","SPL",
                  paste("MFI",analytes,sep="."), paste("CTS",analytes,sep="."), paste("RES",analytes,sep="."))]

    # Convert calibrators and controls names if necessary
    cals    = z.run[z.run$Type %in% c("Standard","QC"),]
    bkg     = vector(mode="numeric", length=length(analytes))
    for (i1 in 1:nrow(cals)){
        if (nchar(cals$SPL[i1])>3 & cals$Type[i1]=="Standard") 
            cals$SPL[i1] = paste("St", substr(cals$SPL[i1], nchar(cals$SPL[i1]), nchar(cals$SPL[i1])+1), sep="")
        if (nchar(cals$SPL[i1])>4 & cals$Type[i1]=="QC") 
            cals$SPL[i1] = paste("Con", substr(cals$SPL[i1], nchar(cals$SPL[i1]), nchar(cals$SPL[i1])+1), sep="")
        }
        
    # Add calibrators & controls
    for (i1 in 1:length(analytes)) {
        for (i2 in c("St1","St2","St3","St4","St5","St6","St7","St8","St9","ConA","ConB")) {
            z.run[z.run$Loc %in% cals$Loc[cals$SPL==i2], paste("conc", analytes[i1], sep=".")] = 
                cl.r[cl.r$Analyte== analytes[i1] ,i2]
            }
        z.run[z.run$Type=="Background", paste("conc", analytes[i1], sep=".")] = 0
        
        # Remove background
        bkg[i1]  = mean(z.run[z.run$Type=="Background", paste("MFI", analytes[i1], sep=".")], na.rm=TRUE)
        z.run[, paste("MFI", analytes[i1], sep=".")] = z.run[, paste("MFI", analytes[i1], sep=".")] - bkg[i1]

        # Protect against counts <20
        z.run[z.run$SPL==i2 & z.run[, paste("CTS", analytes[i1], sep=".")]<20, paste("conc", analytes[i1], sep=".")] = NA
        }

    # Final fixes
    z.run$SPL = factor(z.run$SPL)
    #z.run$Type[z.run$Type=="Background"] = "Missing"
        
    # The end
    return(structure(z.run, file = file,
        Assay      = c(Platform   = as.character(z.info$platform),
                       SN         = as.character(z.info$sn),
                       Session    = as.character(z.info$ses)),
        Lot        = c(Lot        = as.character(z.info$lot), 
                       N          = z.info$n), 
        Analytes   = as.character(analytes),
        Units      = as.character(units),
        calibrators= cl.r,
        coefs      = z.COEF,
        fit        = z.FIT,
        Date       = as.character(as.Date(as.character(z.info$date), format="%m/%d/%Y %H:%M")),
        Operator   = as.character(z.info$operator),
        Background = bkg,
        class      = c("ima","data.frame")))
}

