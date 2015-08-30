#!/usr/bin/env python

import os
import random
import string
import gdal
import math
from itertools import product

import shapefile
import numpy as np
import pandas as pd
import pyproj
import shapely.geometry
import spatial_csv_to_kml


# 1. Pre-defined functions

# 1.1 Creating sampling squares


def create_square(x, y, pixelWidth, pixelHeight):
    """
    Returns a square as an instance of shapely.geometry.Polygon.
    Args:
        x, y: center of the square
        side: side length of the square
    """
    return shapely.geometry.Polygon([(x, y),
                                     (x + pixelWidth, y),
                                     (x + pixelWidth, y + pixelHeight),
                                     (x + pixelWidth, y + pixelHeight)])


# 1.2 Calculating GID

def GID(location):
    res_pixel = 1000
    xgid = int(math.floor(math.fabs(location[0]) / res_pixel))
    ygid = int(math.floor(math.fabs(location[1]) / res_pixel))
    if location[0] < 0:
        gidx = "W" + str(xgid)
    else:
        gidx = "E" + str(xgid)
    if location[1] < 0:
        gidy = "S" + str(ygid)
    else:
        gidy = "N" + str(ygid)

    GID = gidx + "-" + gidy
    return GID


# 2. Create the ROI cropland prediction map

# first, download geosurvey cropland prediction map from
# https://www.dropbox.com/s/209iemv4b2jla9f/geosurvey_h2o_crp_predictions.tif?dl=0
#
#
# then, create a shape file for your ROI projected in '+proj=laea +lat_0=5
# +lon_0=20 +ellps=WGS84 +units=m +no_defs'

# #### 2.1 Using the ROI shapefile (projected in '+proj=laea +lat_0=5 +lon_0=20 +ellps=WGS84 +units=m +no_defs') clipping geosurvey cropland prediction map

shpfile = "roi_shp/bako_8k_laea.shp"
inputfile = "geosurvey_h2o_crp_predictions.tif"
cmd = "gdalwarp -cutline " + shpfile + \
    " -crop_to_cutline -srcnodata \"nan\" -dstnodata \"nan\" " + inputfile + " output.tif"
os.system(cmd)


# #### 2.2 Find pixels with presence probability larger than cutoff (here set 0.5) in geosurvey crop prediction map

cutoff = 0.5

crpdat = gdal.Open("output.tif", gdal.GA_ReadOnly)
crp_prob_np = np.asarray(crpdat.GetRasterBand(1).ReadAsArray())
originX, pixelWidth, rx, originY, ry, pixelHeight = crpdat.GetGeoTransform()

crp_presence_loc = np.where(crp_prob_np > cutoff)
n_presence = crp_presence_loc[0].shape[0]


# #### delete the temporary file "output.tif"


cmd = "rm output.tif"
os.system(cmd)


# #### 2.3 ROI boundary

shp = shapefile.Reader(shpfile, 'rb')
shape = shp.shapes()[0]
polygon = shapely.geometry.asShape(shape)


# ### 3. Start Sampling

# #### 3.1 set number of 1k by 1k grids (n) and number of 100m by 100m grid centers for each 1k by 1k grid (n_pixel) to be sampled

res_1 = 1000
res_2 = 100
n = 10
n_pixel = 10


# #### 3.2 Sampling Start

sampled_locs = []
sampled_locs_idx = []
sampled_1k_locs = []
GID_1k = []
GID_100m = []
i = 0
while i < n:
    init_flag = True
    """
    only keep the squares which are completely contained in the polygon
    """
    while init_flag or not polygon.contains(
            pixel_current) or sampled_loc_idx in sampled_locs_idx:
        init_flag = False
        sampled_loc_idx = random.sample(xrange(n_presence), 1)
        x_current_origin = originX + pixelWidth * \
            crp_presence_loc[0][sampled_loc_idx[0]]
        y_current_origin = originY + pixelHeight * \
            crp_presence_loc[1][sampled_loc_idx[0]]
        pixel_current = create_square(
            x_current_origin, y_current_origin, res_1, res_1)

    sampled_locs_idx.extend(np.repeat(sampled_loc_idx, n_pixel))
    sampled_1k_locs.extend([[x_current_origin, y_current_origin]
                            for x in xrange(n_pixel)])

    GID_1k.extend(
        np.repeat(GID([x_current_origin, y_current_origin]), n_pixel))

    """
    generate 100m by 100m grids
    """
    grid_100m = map(lambda u: [u[0] + x_current_origin,
                               u[1] + y_current_origin],
                    list(product(xrange(-res_1 / 2 + 100 / 2,
                                        res_1 / 2,
                                        100),
                                 xrange(-res_1 / 2 + 100 / 2,
                                        res_1 / 2,
                                        100))))
    sampled_100m_idx = random.sample(xrange(100), n_pixel)
    GID_100m.extend([GID([x_current_origin, y_current_origin]
                         ) + "-" + str(u) for u in sampled_100m_idx])
    sampled = [grid_100m[u] for u in sampled_100m_idx]
    sampled_locs.extend(sampled)
    i += 1


# #### 3.3 project LAEA coordinates back to LL

project_string = '+proj=laea +lat_0=5 +lon_0=20 +ellps=WGS84 +units=m +no_defs'
p = pyproj.Proj(project_string)
sampled_locs_latlon = map(lambda x: list(
    p(x[0], x[1], inverse=True)), sampled_locs)
sampled_1k_locs_latlon = map(lambda x: list(
    p(x[0], x[1], inverse=True)), sampled_1k_locs)


# ### 4. output

output = [tuple(GID_1k)]
output.extend(zip(*sampled_locs_latlon))
output.extend([tuple(GID_100m)])
output.extend(zip(*sampled_1k_locs_latlon))
output_pd = pd.DataFrame(output)
output_pd = output_pd.transpose()
output_pd.columns = ['GID_1k', 'x', 'y', 'GID_100m', 'x_1k', "y_1k"]

if not os.path.exists('output'):
    os.makedirs('output')
if not os.path.exists('output/drone_flight_1k'):
    os.makedirs('output/drone_flight_1k')
output_pd.to_csv('output/sampled_loc.csv', index=False)

spatial_csv_to_kml.csv_to_kml('output/sampled_loc.csv')
