
# coding: utf-8

# # Ghana Soil Information Service Sampling Plan (2015)
# 
# Jiehua Chen (jc3288@columbia.edu)

# In[3]:

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


# ### Getting Geosurvey data: those data are shared through dropbox link: 
# geosurvey data: https://www.dropbox.com/s/f0fyy3dlh44ar0b/geosurvey_results.zip?dl=0
# 
# shapefile: https://www.dropbox.com/s/fj61a0rvxyji23d/shapefile.zip?dl=0

# In[19]:

shpfile = "data/ROI_GhaSIS_GAM_districts/ROI_GhaSIS_GAM_laea.shp"
inputfile = "data/geosurvey_crp_prediction_10k.tif"
cmd = "gdalwarp -cutline " + shpfile + " -crop_to_cutline -srcnodata \"nan\" -dstnodata \"nan\" " + inputfile + " mean_output.tif"
os.system(cmd)

inputfile = "data/geosurvey_h2o_crp_predictions.tif"
cmd = "gdalwarp -cutline " + shpfile + " -crop_to_cutline -srcnodata \"nan\" -dstnodata \"nan\" " + inputfile + " mean_output_1k.tif"
os.system(cmd)


# In[20]:

crpdat = gdal.Open("mean_output.tif", gdal.GA_ReadOnly)
crp_prob_np = np.asarray(crpdat.GetRasterBand(1).ReadAsArray())
originX, pixelWidth, rx, originY, ry, pixelHeight = crpdat.GetGeoTransform()

crpdat_1k = gdal.Open("mean_output_1k.tif", gdal.GA_ReadOnly)
crp_prob_np_1k = np.asarray(crpdat.GetRasterBand(1).ReadAsArray())
originX_1k, pixelWidth_1k, rx_1k, originY_1k, ry_1k, pixelHeight_1k = crpdat.GetGeoTransform()


# In[7]:

cutoff_mean = 0.31
crp_presence_loc = np.where((crp_prob_np>cutoff_mean))
n_presence = crp_presence_loc[0].shape[0]
n_presence


# In[8]:

def getcoords(idxx, idxy, w=pixelWidth, h=pixelHeight, x0 = originX, y0=originY):
    x = x0 + w * idxy
    y = y0 + h * idxx
    return [x,y]

crp_presence_coords = list(starmap(getcoords, zip(*crp_presence_loc)))


# In[23]:

shpfile = "data/ROI_GhaSIS_GAM_districts/ROI_GhaSIS_GAM_laea.shp"
districts_roi_shp = shapefile.Reader(shpfile, 'rb')
districts_roi_names = map(lambda x: x[6], districts_roi_shp.records())


# ### Find locations within each district, which have cropland presence probability larger than the cutoff value

# In[24]:

regions = districts_roi_shp.shapes()
square_size = 10000
progress_mark = 500


# In[25]:

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


# In[26]:

points_with_regions_pd = pd.DataFrame(zip(*points_with_regions)).transpose()
points_with_regions_pd.columns = ['x', 'y','district_idx']
points_with_regions_pd.to_csv("locations_crop_10k.csv", index=False)


# ### Find grids with high cropland presence and their locations

# In[27]:

points_with_regions_pd = pd.DataFrame(zip(*points_with_regions)).transpose()
points_with_regions_pd.columns = ['x','y','district']
district_highcrp = map(lambda x: int(x), points_with_regions_pd.groupby('district').count().index)
district_highcrp_count = points_with_regions_pd.groupby('district').count()
crp_area_prob = [district_highcrp_count[district_highcrp_count.index==k].x.get_values()[0] if k in district_highcrp else 0 for k in xrange(len(districts_roi_names))]
crp_area_prob = map(lambda x: x*1.0/sum(crp_area_prob), crp_area_prob)


# ### sample size

# In[36]:

n_10k = 100
n_1k = 10
n_100m = 5
n_10k_perdistrict = map(lambda x: int(math.ceil(n_10k*x)), crp_area_prob)
sum(n_10k_perdistrict)


# In[39]:

random.seed(20150906)
reload(multistage_sampling)


# ### Start Sampling: the results are sampling results of csv files for each district, csv, kml files for each 10k-by-10k, kml files for drone flights organized in the corresponding district folder

# In[41]:

for k, sample_n_10k in enumerate(n_10k_perdistrict):
    if sample_n_10k >0:
        district_locs = points_with_regions_pd[points_with_regions_pd.district==k]
        sampled_locs_idx = random.sample(list(district_locs.index), sample_n_10k)
        sampled_locs =  district_locs.ix[sampled_locs_idx]
        x_current_lower = list(sampled_locs.x.get_values()+5000)
        y_current_lower = list(sampled_locs.y.get_values()-5000)
        multistage_sampling.sample(x_current_lower,  y_current_lower, n_100m, districts_roi_names[k], "mean_output_1k.tif", 0.51, number_of_1_by_1=n_1k)
    


# In[49]:

import glob
import csv
current_dir = "/Users/jiehuachen/Documents/research/afsis/sampling_protocol/sampling/Ghana/output/"


# In[54]:

os.chdir(current_dir)
district_folders = glob.glob('*')
os.chdir(os.path.join(current_dir, district_folders[0])) 
csvfile=glob.glob("*.csv")[0]
locs_data = pd.read_csv(csvfile)
locs_data.keys().get_values()
header = locs_data.keys().get_values()
os.chdir(current_dir)
with open(current_dir+'sampled_locs_total.csv', 'w+') as csvfile:
    headerwriter = csv.writer(csvfile, delimiter=',')
    headerwriter.writerow(header)
csvfile.close()

with open(current_dir+"sampled_locs_total.csv", 'a') as f: 
    for folder in district_folders:
        os.chdir(os.path.join(current_dir, folder)) 
        for csvfile in glob.glob("*.csv"):
            locs_data = pd.read_csv(csvfile)
            locs_data['district_name'] = folder
            locs_data.to_csv(f, header=False, index=False)

