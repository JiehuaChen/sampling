
source("mstage.R")
mstage_sampling <- function(location, n_pixel){

    utm_zone = ceiling(1 + floor((location[1]+180)/6.0))
    xycenter.location <- project(as.matrix(location),paste("+proj=utm +zone=", utm_zone,  " +datum=WGS84 +units=m +no_defs", sep=""))
    xoff <-  xycenter.location[1]-5000
    yoff <- xycenter.location[2]-5000

    ### xdim & ydim specifiy the number of cells to be sampled in x & y directions
    xdim <- 100
    ydim <- 100

    ### specify grid resolution (grain, in m)
    grain <- 100

    ### generate the grid
    grid <- as.data.frame(coordinates(GridTopology(c(xoff,yoff), c(grain,grain), c(xdim,ydim))))
    colnames(grid) <- c("x", "y")

    ### set up level ID's at the desired scales (res.pixel, in m)
    res.pixel <- c(5000, 1000, grain) 
    X2 <- ceiling((grid$x-xoff+1/2*xdim)/res.pixel[1])
    Y2 <- ceiling((grid$y-yoff+1/2*ydim)/res.pixel[1])
    L2 <- cleanstrata(paste(X2, Y2, sep=""))

    X1 <- ceiling((X2*res.pixel[1]-(grid$x-xoff))/res.pixel[2])
    Y1 <- ceiling((Y2*res.pixel[1]-(grid$y-yoff))/res.pixel[2])
    L1 <- cleanstrata(paste(X1, Y1, sep=""))

    X0 <- ceiling((X1*res.pixel[2]+(grid$x-xoff-X2*res.pixel[1]))/res.pixel[3])
    Y0 <- ceiling((Y1*res.pixel[2]+(grid$y-yoff-Y2*res.pixel[1]))/res.pixel[3])
    L0 <- cleanstrata(paste(X0, Y0, sep=""))

    ### update grid with level ID's
    grid <- cbind(grid, L2, L1, L0)

    ### draw a sample
    ### s.size specifies the sample size at each level (e.g. L2=4, L1=4, L0=10)
    s.size <- list(4, rep(4, 4), rep(n_pixel, 16))
    s <- mstage(grid, stage=list("cluster", "cluster", ""), varnames=list("L2", "L1", "L0"), size=s.size, method="srswor", description=TRUE)

    sample.DGG <- getdata(grid, s[[3]])
    #### plot1 <- xyplot(y~x, data=sample.DGG, xlab="Easting (m)", ylab="Northing (m)", pch=3, cex=0.5, asp=1)

    ### define a spatial point data object & project from DGG to geographic coordinates 
    DGG <- CRS(paste("+proj=utm +zone=",utm_zone, " +datum=WGS84 +units=m +no_defs", sep=""))

    coordinates(sample.DGG) <- ~x+y
    proj4string(sample.DGG) <- DGG
    sample.LL <- spTransform(sample.DGG, CRS("+proj=longlat +datum=WGS84"))

    ### write an OGR file (e.g. KML) for visualization, navigation ... 
    ### for different OGR drivers see: http://www.gdal.org/ogr/ogr_formats.html
    write.csv(as.data.frame(sample.LL), paste("csv/E", location[1],"N", location[2], "_", n_pixel*16, ".csv", sep=""))
}
