.required <- c("graphics", "gplots", "methods", "stats", "plotrix")

cv <- function(x) {
    mu = mean(x, na.rm=TRUE)
    sg = sd(x, na.rm=TRUE)
    return(sg/mu*100)
}

mesdci <- function(a, dig=4, shown=FALSE) {
    n       = length(na.omit(a))
    mean2   = format(mean(a, na.rm=TRUE), digits=dig)
    sd2     = format(sd(a, na.rm=TRUE), digits=dig)
    cv2     = round(sd(a, na.rm=TRUE)/mean(a, na.rm=TRUE)*100, 2)
    median2 = format(median(a, na.rm=TRUE), digits=dig)
    lci     = format(quantile(a, probs=0.025, na.rm=TRUE), digits=dig)
    uci     = format(quantile(a, probs=0.975, na.rm=TRUE), digits=dig)
    if (shown) {
        all     = c(n,mean2,sd2,cv2,median2,lci,uci)
        names(all) = c("N","Mean","SD","CV","Median","LCI","UCI") }
    else {
        all     = c(mean2,sd2,cv2,median2,lci,uci)
        names(all) = c("Mean","SD","CV","Median","LCI","UCI") }
    return(all)
}
 