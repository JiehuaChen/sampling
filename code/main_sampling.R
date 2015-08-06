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

sample_locs <- mstage_sampling(center_loc, 10, 8000)

# generate kml files
system("python  spatial_csv_to_kml.py")

# generate gpx file for GARMIN navigation
wpts <- sample_locs[, c("x", "y", "GID_100m")]

write.csv(wpts, "csv/Waypoints.csv", row.names=FALSE)
# note that this system call has to point to the location of your GPSBabel application
system("gpsbabel -i csv -f csv/Waypoints.csv -o gpx -F csv/Waypoints.gpx")

