
## Multistage Sampling Based on Geosurvey Results

#### Jiehua Chen, Markus Walsh, William Wu

### Load libraries


library(lattice)
require(proj4)
require(sampling)
library(MASS)
library(rgdal)
library(raster)


### load multistage sampling code


source("../code/mstage_sampling.R")


### create output directory


dir.create("csv")
fileincsv <- list.files("csv",full.names=TRUE )
for(i in fileincsv){
    print(i)
    file.remove(i)
}


### Geosurvey Data
#### csv file here is shared in: https://www.dropbox.com/s/n6gyeohpktkx060/1MQ_CRP_pos.csv?dl=0
#### tif file is shared in: https://www.dropbox.com/s/209iemv4b2jla9f/geosurvey_h2o_crp_predictions.tif?dl=0

geosurvey_data_file = "1MQ_CRP_pos.csv"
crp_map = "geosurvey_h2o_crp_predictions.tif"


### sampling procedure


#### geosurvey_data_file: file name of geosurvey data
#### crp_map: crpland probability map
#### the cutoff value is the threshold for cropland presence
#### colname, response: we only randomly choose from rows with the values of column with column name being colname is "response"
#### for example: colname = 'Cropland.present'; response = 'Yes' meaning that we select among cropland pixels
#### n_pixel: number of points selected within each 1k by 1k grid
#### latcolname and loncolname are the names for columns with latitude and longitude values
#### pixelWith and pixelHeight is the pixel size for geosurvey
#### size is the length of UUID


locsfromcsv <- function(geosurvey_data_file, crp_map, cutoff, n, n_pixel, colname, response, latcolname, loncolname){
    row_idx = 0
    geosurvey_data = read.csv(geosurvey_data_file)
    crpland_data = geosurvey_data[geosurvey_data[colname]==response, ]
    
    ##### load in cropland presense map
    crpmap <- stack(crp_map)
    k <- 0
    while (k<n){
        ##### each data point is sampled with each probability 
        sampled_loc_idx = sample(1:dim(crpland_data)[1], 1)
        sampled_center = crpland_data[sampled_loc_idx,][, c(loncolname, latcolname)]
        project_string = '+proj=laea +lat_0=5 +lon_0=20 +ellps=WGS84 +units=m +no_defs' 
        loc = project(as.matrix(sampled_center), project_string) 
        x_origin <- loc[1]-5000
        y_origin <- loc[2]+5000

        x_pixels <- seq(x_origin, by=1000, length=10)
        y_pixels <- seq(y_origin, by=-1000, length=10)
        xy_pixels <- cbind(rep(x_pixels, 10), rep(y_pixels, each=10))
        crpdata_grid <- extract(crpmap, xy_pixels)
        #### check if the whole 10K by 1K grid centered by the selected location have high probability of cropland presence
        if(mean(crpdata_grid>cutoff, na.rm=TRUE)==1){
            #### the multistage sampling can incorporate accessibility measurement later
            mstage_sampling(sampled_center, n_pixel)
            k <- k+1
            print(k)
        }
    }
}



colname="CRP"
response ="Y"
loncolname = "Lon";
latcolname = "Lat";



locsfromcsv(geosurvey_data_file, crp_map, 0.26, 3, 10, colname, response, latcolname, loncolname)


### Create KML files


system("python  ../code/spatial_csv_to_kml.py")





