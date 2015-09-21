
# coding: utf-8

# In[1]:

import os
import random
import string
import gdal
import math
from itertools import product
from itertools import starmap
import shapefile
import numpy as np
import pandas as pd
import pyproj
import shapely.geometry
import spatial_csv_to_kml
from rasterstats import zonal_stats

from rtree import index

from shapely.geometry import Polygon
from shapely.geometry import Point

import multistage_sampling 


# ### Getting Geosurvey data: those data are shared through dropbox link: https://www.dropbox.com/s/zrq0bx83dgzpvsp/data.zip?dl=0
# 

# In[2]:

shpfile = "data/tansis_sampling_districts_laea.shp"
inputfile = "data/TZ_mean_10k.tif"
cmd = "gdalwarp -cutline " + shpfile + " -crop_to_cutline -srcnodata \"nan\" -dstnodata \"nan\" " + inputfile + " mean_output.tif"
os.system(cmd)

inputfile = "data/TZ_crp_ens_1k.tif"
cmd = "gdalwarp -cutline " + shpfile + " -crop_to_cutline -srcnodata \"nan\" -dstnodata \"nan\" " + inputfile + " mean_output_1k.tif"
os.system(cmd)


# In[4]:

crpdat = gdal.Open("mean_output.tif", gdal.GA_ReadOnly)
crp_prob_np = np.asarray(crpdat.GetRasterBand(1).ReadAsArray())
originX, pixelWidth, rx, originY, ry, pixelHeight = crpdat.GetGeoTransform()

crpdat_1k = gdal.Open("mean_output_1k.tif", gdal.GA_ReadOnly)
crp_prob_np_1k = np.asarray(crpdat.GetRasterBand(1).ReadAsArray())
originX_1k, pixelWidth_1k, rx_1k, originY_1k, ry_1k, pixelHeight_1k = crpdat.GetGeoTransform()


# In[5]:

cutoff_mean = 0.7
cutoff_sd = 0.1
crp_presence_loc = np.where((crp_prob_np>cutoff_mean))
n_presence = crp_presence_loc[0].shape[0]
n_presence


# In[6]:

def getcoords(idxx, idxy, w=pixelWidth, h=pixelHeight, x0 = originX, y0=originY):
    x = x0 + w * idxy
    y = y0 + h * idxx
    return [x,y]

crp_presence_coords = list(starmap(getcoords, zip(*crp_presence_loc)))


# In[7]:

shpfile = "data/tansis_sampling_districts_laea.shp"
districts_roi_shp = shapefile.Reader(shpfile, 'rb')
districts_roi_names = map(lambda x: x[6], districts_roi_shp.records())


# ### Find locations within each district, which have cropland presence probability larger than the cutoff value

# In[8]:

regions = districts_roi_shp.shapes()
square_size = 10000
progress_mark = 500


# In[9]:

points_with_regions = []

rtree_idx = index.Index()

for i,r in enumerate(regions):
    rtree_idx.insert(i,r.bbox)
    
n = len(crp_presence_coords)
locs_inregion = []
for k,georef in enumerate(crp_presence_coords):
    if 0 == k % progress_mark: print "{}/{} ~= {:.2f}%% complete".format(k,n,100*float(k)/n)
    x,y = georef
    square = shapely.geometry.geo.box(x, y-square_size, x+square_size, y)
    for j in rtree_idx.intersection((x,y-square_size,x+square_size,y)):
        if shapely.geometry.asShape(regions[j]).contains(square):
            points_with_regions.append([x,y,j]) 
            break # WARNING: this assumes exactly one region will contain the point


# In[10]:

points_with_regions_pd = pd.DataFrame(zip(*points_with_regions)).transpose()
points_with_regions_pd.columns = ['x', 'y','district_idx']
points_with_regions_pd.to_csv("locations_crop_10k.csv", index=False)


# ### Find high cropland presence grids and their locations

# In[11]:

points_with_regions_pd = pd.DataFrame(zip(*points_with_regions)).transpose()
points_with_regions_pd.columns = ['x','y','district']
district_highcrp = map(lambda x: int(x), points_with_regions_pd.groupby('district').count().index)
district_highcrp_count = points_with_regions_pd.groupby('district').count()
crp_area_prob = [district_highcrp_count[district_highcrp_count.index==k].x.get_values()[0] if k in district_highcrp else 0 for k in xrange(len(districts_roi_names))]
crp_area_prob = map(lambda x: x*1.0/sum(crp_area_prob), crp_area_prob)


# ### sample size

# In[12]:

n_10k = 200
n_1k = 4
n_100m = 5
n_10k_perdistrict = map(lambda x: int(math.ceil(n_10k*x)), crp_area_prob)
sum(n_10k_perdistrict)


# ### Start Sampling: the results are sampling results of csv files for each district, csv, kml files for each 10k-by-10k, kml files for drone flights organized in the corresponding district folder

# In[28]:

random.seed(20150906)
reload(multistage_sampling)


# In[29]:

for k, sample_n_10k in enumerate(n_10k_perdistrict):
    if sample_n_10k >0:
        district_locs = points_with_regions_pd[points_with_regions_pd.district==k]
        sampled_locs_idx = random.sample(list(district_locs.index), sample_n_10k)
        sampled_locs =  district_locs.ix[sampled_locs_idx]
        x_current_lower = list(sampled_locs.x.get_values()+5000)
        y_current_lower = list(sampled_locs.y.get_values()-5000)
        multistage_sampling.sample(x_current_lower,  y_current_lower, n_100m, districts_roi_names[k], "mean_output_1k.tif")
    
        sampled_data = pd.read_csv('output/'+districts_roi_names[k]+'/'+districts_roi_names[k]+'_'+str(0)+'.csv')[['y','x']]
        total_file = 'output/'+districts_roi_names[k]+'/'+districts_roi_names[k]+'.csv'
        sampled_data.to_csv(total_file, index=False)
        if sample_n_10k>1:
            with open(total_file, 'a') as f:            
                for s in xrange(1, sample_n_10k):
                    filename = 'output/'+districts_roi_names[k]+'/'+districts_roi_names[k]+'_'+str(s)+'.csv'
                    sampled_data = pd.read_csv(filename)[['y','x']]
                    sampled_data.to_csv(f, header=False, index=False)
            
        filenames_gpx = [os.system('rm '+'"output/'+districts_roi_names[k]+'/'+districts_roi_names[k]+'_'+str(s)+'_Waypoints'+'.csv"') for s in xrange(sample_n_10k)]
        


# In[30]:

os.chdir("output/")
import glob
locs_csvfiles =  glob.glob("*")
locs_data = pd.read_csv(os.path.join(locs_csvfiles[0],locs_csvfiles[0]+".csv"))
locs_data['district_name'] = locs_csvfiles[0]
locs_data.to_csv("sampled_locs_total.csv", index=False)
with open("sampled_locs_total.csv", 'a') as f: 
    for locs_csvfile in locs_csvfiles[1:]:
        sampled_data = pd.read_csv(os.path.join(locs_csvfile, locs_csvfile+".csv"))
        sampled_data['district_name'] = locs_csvfile
        tmp = [locs_data, sampled_data]
        locs_data = pd.concat(tmp)
        sampled_data.to_csv(f, header=False, index=False)

