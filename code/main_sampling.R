
source("mstage_sampling.R")

### create output directory


dir.create("csv")
fileincsv <- list.files("csv",full.names=TRUE )
for(i in fileincsv){
    print(i)
    file.remove(i)
}

center_loc <- matrix(c(36, -3), 1, 2)

system("python  spatial_csv_to_kml.py")

