
# coding: utf-8

# In[73]:

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

#### data folder is zipped and shared in https://www.dropbox.com/s/zrq0bx83dgzpvsp/data.zip?dl=0

shpfile = "data/sagcot_districts_laea.shp"
inputfile = "data/geosurvey_crp_prediction_10k.tif"
cmd = "gdalwarp -cutline " + shpfile + " -crop_to_cutline -srcnodata \"nan\" -dstnodata \"nan\" " + inputfile + " output.tif"
os.system(cmd)



crpdat = gdal.Open("output.tif", gdal.GA_ReadOnly)
crp_prob_np = np.asarray(crpdat.GetRasterBand(1).ReadAsArray())
originX, pixelWidth, rx, originY, ry, pixelHeight = crpdat.GetGeoTransform()



shpfile = "sagcot_districts_laea.shp"
districts_roi_shp = shapefile.Reader(shpfile, 'rb')
districts_roi_names = map(lambda x: x[6], districts_roi_shp.records())


### Find locations within each district, which have cropland presence probability larger than the cutoff value


def sum_narm(x, cutoff=0.5):
    mdat = np.ma.masked_array(x>cutoff,np.isnan(x))
    mm = np.sum(mdat)
    return mm
def highlocs(x, cutoff=0.5):
    mdat2 = np.ma.masked_outside(x, cutoff, 1)
    return np.where(mdat2)
stats = zonal_stats("sagcot_districts_laea.shp", "output.tif", add_stats={'sum_narm':sum_narm, 'highlocs':highlocs})
locs_crp = map(lambda x: [x['highlocs'][0].compressed(), x['highlocs'][1].compressed()], stats)
crp_area = map(lambda x: x['sum_narm'], stats)
crp_area_prob = map(lambda x: x/ float(sum(crp_area)), crp_area)



cmd = "rm output.tif"
os.system(cmd)


# ### sample size

n_10k = 100
n_1k = 4
n_100m = 5
n_10k_perdistrict = map(lambda x: int(math.ceil(n_10k*x)), crp_area_prob)


# ### Start Sampling


for k, sample_n_10k in enumerate(n_10k_perdistrict):
    sampled_locs_idx = random.sample(xrange(crp_area[k]), sample_n_10k)
    x_current_lower = [originX + pixelWidth *locs_crp[k][0][i]- pixelWidth/2 for i in  sampled_locs_idx]
    y_current_lower = [originY + pixelHeight *locs_crp[k][1][i] - pixelHeight/2 for i in  sampled_locs_idx]
    multistage_sampling.sample(x_current_lower,  y_current_lower, n_100m, districts_roi_names[k])


