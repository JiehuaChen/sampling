   
GID <- function(location){
    res.pixel <- 1000
    xgid <- floor(location[,1]/res.pixel)
    ygid <- floor(location[,2]/res.pixel)
    gidx <- ifelse(location[,1]<0, paste("W", xgid, sep=""), paste("E", xgid, sep=""))
    gidy <- ifelse(location[,2]<0, paste("S", ygid, sep=""), paste("N", ygid, sep=""))
    GID <- paste(gidx, gidy, sep="-")
    return(list(xgid=xgid, ygid=ygid, GID=GID))
}
