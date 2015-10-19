
# coding: utf-8

# In[1]:

library(rgdal)


# In[65]:

geosurvey_answers <- read.csv("export-tansis_otheregions.csv", stringsAsFactors=FALSE)
sampled_locs <- read.csv("sampled_locs_total.csv", stringsAsFactors=FALSE)


# In[66]:

geosurvey_locs_laea <- project(cbind(geosurvey_answers$Longitude, geosurvey_answers$Latitude), '+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs')
locs_laea <- project(cbind(sampled_locs$x, sampled_locs$y), '+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs')
locinfo <- matrix(NA, 1, dim(sampled_locs)[2])
colnames(locinfo) <- names(sampled_locs)


# In[74]:

for(i in 1:dim(geosurvey_answers)[1]){
    dist_tmp <- sqrt((geosurvey_locs_laea[i, 1]-locs_laea[,1])^2 + (geosurvey_locs_laea[i, 2]-locs_laea[,2])^2)
    locinfo <- rbind(locinfo, sampled_locs[which.min(dist_tmp), ])
}


# In[75]:

locinfo <- locinfo[-1, ]


# In[76]:

geosurvey_answers <- data.frame(geosurvey_answers, locinfo, x=geosurvey_locs_laea[,1], y=geosurvey_locs_laea[,2])


# In[77]:

dim(locinfo)


# In[78]:

table(geosurvey_answers_yes$GID_1k)


# In[79]:

geosurvey_answers_yes <- geosurvey_answers[geosurvey_answers$Cropland.Present.=="Yes",]


# In[80]:

str(geosurvey_answers_yes)


# In[81]:

sampled_random_locs <- geosurvey_answers_yes[sample(1:dim(geosurvey_answers_yes)[1], dim(geosurvey_answers)[1]/2), ]


# In[58]:

sampled_random_locs <- sampled_random_locs[,c("x", "y", "GID_100m", "x_1k", "y_1k","GID_1k", "district")]


# In[59]:

str(sampled_random_locs)


# In[60]:

district_names <- unique(sampled_random_locs$district)
district_names


# In[61]:

for(district in district_names){
    dir.create(file.path("output_phaseII/", district))
    write.csv(sampled_random_locs[sampled_random_locs$district==district, ], file.path("output_phaseII/", district, paste(district, ".csv", sep="")), row.names=FALSE)       
}




