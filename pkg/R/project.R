  
project.means      <- function(x) {
    minus = match(c("ID","Type","Loc","SPL","Lot"), names(x))
    res  = aggregate(x[, -minus], by=list(ID=x$ID,Type=x$Type,SPL=x$SPL, Lot=x$Lot), FUN=mean, na.rm=TRUE)
    res  = res[res$Type!="Missing",]
    res  = res[order(res$ID, res$Type, res$SPL),]
    
    res$RID	    = res$VISIT = NA
    res$RID[res$Type=="Sample"]   = as.numeric(substr(as.character(res$SPL[res$Type=="Sample"]),0,4))
    res$VISIT[res$Type=="Sample"] = substr(as.character(res$SPL[res$Type=="Sample"]),5,8)
    res$VISIT[!is.na(res$VISIT) & res$VISIT=="v4-n"] = "v4"
    res$VISIT	= ordered(res$VISIT, levels=c("v2","v4","v6","v7","v8"), labels=c("BL","M12","M24","M36","M48"))
    
    invisible(res)
}

project.qcplot     <- function(x, samples=c("ConA","ConB"), type="pred", breaks=c(1,length(unique(x$ID))), log=FALSE, ...) {
    z.analytes = substr(names(x)[substr(names(x),0,4)=="MFI."],5,50)
    n  		   = length(z.analytes)
    for (i3 in samples) {
        my  = max(x[x$SPL==i3, paste(rep(type, n), z.analytes, sep=".")], na.rm=TRUE)
        mmy = min(x[x$SPL==i3, paste(rep(type, n), z.analytes, sep=".")], na.rm=TRUE)
        mx  = nrow(x[x$SPL==i3,])
        plot(0,0, ylim=c(ifelse(log==TRUE,mmy*0.9,0), my), xlim=c(1, mx), xlab = "Plate #", 
            ylab= paste(ifelse(log==TRUE, "log(",""), ifelse(type=="pred","Concentration [pg/mL]",type), ifelse(log==TRUE, ")",""), sep=""),
            type="n", main= paste("Stability -", i3), log=ifelse(log==TRUE,"y",""))

        for (i1 in 1:length(z.analytes)) {
            i2 = paste(type, z.analytes[i1], sep=".")
            for (i4 in 1:(length(breaks)-1)) {
                a.r = c(breaks[i4], breaks[i4+1]) - c(0, ifelse(i4<length(breaks)-1,1,0))
                a.m = mean(x[x$SPL==i3, i2][a.r[1]:a.r[2]], na.rm=TRUE)
                a.s = sd(x[x$SPL==i3, i2][a.r[1]:a.r[2]], na.rm=TRUE)
                if (is.na(a.s)) a.s = 0
                a.r = a.r + c(-0.5, 0.5)
                for(i5 in 1:2) {
                    polygon(rbind(cbind(c(a.r[1], a.r[2]), rep(a.m-i5*a.s,2)),
                        cbind(c(a.r[2], a.r[1]), rep(a.m+i5*a.s,2))),
                        col=c("#ff873320","#8733ff20","#33ff8720")[i1], border=NA)
                    }
                lines(cbind(a.r, rep(a.m,2)), lwd=2, col=c("#ff8733cc","#8733ffcc","#33ff87cc")[i1])
                lines(cbind(a.r, rep(a.m-2*a.s,2)), lwd=0.2, col=c("#ff8733cc","#8733ffcc","#33ff87cc")[i1])
                lines(cbind(a.r, rep(a.m+2*a.s,2)), lwd=0.2, col=c("#ff8733cc","#8733ffcc","#33ff87cc")[i1])
                a.t = as.numeric(mesdci(x[x$SPL==i3, i2][a.r[1]:a.r[2]], dig=3)[c(1,3)])
                pos = ifelse(a.m/my<0.5,-1,1)
                if ((4*a.s)/my>0.33) pos=pos/1.5
                text(a.r[1]+0.5, a.m - 3*a.s*pos, adj=c(0,ifelse(pos==1,1,-0.3)), 
                    labels=paste(z.analytes[i1], paste(" Mean:", formatC(a.t[1], digits=3, format="fg")), 
                    paste( " %CV:  ", round(a.t[2],1),"%", sep=""), sep="\n"), cex=0.8)
                points(x[x$SPL==i3, i2], pch=21, bg= c("#ff873377","#8733ff77","#33ff8777")[i1],
                    col=c("#ff8733","#8733ff","#33ff87")[i1], cex=1)
                }
            }
        legend("bottomleft", z.analytes, fill=c("#ff8733cc","#8733ffcc","#33ff87cc"), bty="n", ncol=3)
        }
    }
