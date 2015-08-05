library(lattice)
require(proj4)
require(sampling)
library(MASS)
library(rgdal)
library(raster)


source("mstage_sampling.R")

### create output directory


dir.create("csv")
fileincsv <- list.files("csv",full.names=TRUE )
for(i in fileincsv){
    print(i)
    file.remove(i)
}

center_loc <- matrix(c(37.115738, 7.042674), 1, 2)

mstage_sampling(center_loc, 10, 8000)

system("python  spatial_csv_to_kml.py")

