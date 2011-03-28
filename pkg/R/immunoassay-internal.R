.required <- c("graphics", "gplots", "methods", "stats")

cv <- function(x) {
    mu = mean(x, na.rm=TRUE)
    sg = sd(x, na.rm=TRUE)
    return(sg/mu*100)
}

mesdci <- function(a, dig=4, shown=FALSE) {
    n       = length(na.omit(a))
    mean    = format(mean(a, na.rm=T), digits=dig)
    sd      = format(sd(a, na.rm=T), digits=dig)
    cv      = round(sd(a, na.rm=T)/mean(a, na.rm=T)*100, 2)
    median  = format(median(a, na.rm=T), digits=dig)
    lci     = format(quantile(a, probs=0.025, na.rm=T), digits=dig)
    uci     = format(quantile(a, probs=0.975, na.rm=T), digits=dig)
    if (shown) {
        all     = c(n,mean,sd,cv,median,lci,uci)
        names(all) = c("N","Mean","SD","CV","Median","LCI","UCI") }
    else {
        all     = c(mean,sd,cv,median,lci,uci)
        names(all) = c("Mean","SD","CV","Median","LCI","UCI") }
    return(all)
}
 