
# coding: utf-8

# In[78]:

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


### Getting Geosurvey data


shpfile = "sagcot_districts_laea.shp"
inputfile = "../data/geosurvey_crp_prediction_10k.tif"
cmd = "gdalwarp -cutline " + shpfile + " -crop_to_cutline -srcnodata \"nan\" -dstnodata \"nan\" " + inputfile + " output.tif"
os.system(cmd)



crpdat = gdal.Open("output.tif", gdal.GA_ReadOnly)
crp_prob_np = np.asarray(crpdat.GetRasterBand(1).ReadAsArray())
originX, pixelWidth, rx, originY, ry, pixelHeight = crpdat.GetGeoTransform()



cutoff = 0.5
crp_presence_loc = np.where(crp_prob_np > cutoff)
n_presence = crp_presence_loc[0].shape[0]



def getcoords(idxx, idxy, w=pixelWidth, h=pixelHeight, x0 = originX, y0=originY):
    x = x0 + w * idxy
    y = y0 + h * idxx
    return [x,y]

crp_presence_coords = list(starmap(getcoords, zip(*crp_presence_loc)))



shpfile = "sagcot_districts_laea.shp"
districts_roi_shp = shapefile.Reader(shpfile, 'rb')
districts_roi_names = map(lambda x: x[6], districts_roi_shp.records())


### Find locations within each district, which have cropland presence probability larger than the cutoff value


regions = districts_roi_shp.shapes()
square_size = 10000
progress_mark = 500



points_with_regions = []

rtree_idx = index.Index()
for i,r in enumerate(regions):
    rtree_idx.insert(i,r.bbox)
    
n = len(crp_presence_coords)
locs_inregion = []
for k,georef in enumerate(crp_presence_coords):
    if 0 == k % progress_mark: print "{}/{} ~= {:.2f}%% complete".format(k,n,100*float(k)/n)
    x,y = georef
    square = shapely.geometry.geo.box(x, y, x+square_size, y+square_size)
    for j in rtree_idx.intersection((x,y,x+square_size,y+square_size)):
        if shapely.geometry.asShape(regions[j]).contains(square):
            points_with_regions.append([x,y,j]) 
            break # WARNING: this assumes exactly one region will contain the point



cmd = "rm output.tif"
os.system(cmd)


### Summarize the high cropland presence locations by districts


points_with_regions_pd = pd.DataFrame(zip(*points_with_regions)).transpose()
points_with_regions_pd.columns = ['x','y','district']
district_highcrp = map(lambda x: int(x), points_with_regions_pd.groupby('district').count().index)
district_highcrp_count = points_with_regions_pd.groupby('district').count()
crp_area_prob = [district_highcrp_count[district_highcrp_count.index==k].x.get_values()[0] if k in district_highcrp else 0 for k in xrange(len(districts_roi_names))]
crp_area_prob = map(lambda x: x*1.0/sum(crp_area_prob), crp_area_prob)


### sample size


n_10k = 100
n_1k = 4
n_100m = 5
n_10k_perdistrict = map(lambda x: int(math.ceil(n_10k*x)), crp_area_prob)


### Start Sampling


for k, sample_n_10k in enumerate(n_10k_perdistrict):
    if sample_n_10k >0:
        district_locs = points_with_regions_pd[points_with_regions_pd.district==k]
        sampled_locs_idx = random.sample(list(district_locs.index), sample_n_10k)
        sampled_locs =  district_locs.ix[sampled_locs_idx]
        x_current_lower = list(sampled_locs.x.get_values())
        y_current_lower = list(sampled_locs.y.get_values())
        multistage_sampling.sample(x_current_lower,  y_current_lower, n_100m, districts_roi_names[k])






