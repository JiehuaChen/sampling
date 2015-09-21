import pyproj
import sys
import math
import numpy as np
import pandas as pd
import random
import spatial_csv_to_kml
import os
import subprocess
import gdal
from itertools import starmap

def cleanstrata(vector):
    """
    Returns list of positions in a sorted vector without duplicates.
    Args:
        vector
    Example:
    > cleanstrata([4, 3, 3]) # > [1, 1, 2]
    """
    a = list(sorted(set(vector)))
    b = {}
    for i in xrange(0, len(a)):
        b[str(a[i])] = i + 1
    return [b[str(v)] for v in vector]


def concatenate(*args, **kwargs):
    """
    Concatenate lists element-wise turning them into strings.
    Args:
        l1, l2, l2.. - arbitrary number of lists to concatenate
        sep - separator in strings
    Returns:
        one list
    Example:
    > concatenate([4, 3, 5], [1, 2, 1], sep=',') # > ['4,1', '3,2', '5,1']
    """
    sep = kwargs.get('sep', '')
    concatenated = []
    for i in xrange(len(args[0])):
        concatenated.append(sep.join([str(a[i]) for a in args]))
    return concatenated

def GID_1k_func(location):
    res_pixel = 1000
    xgid = int(math.floor(math.fabs(location[0])/res_pixel))
    ygid = int(math.floor(math.fabs(location[1])/res_pixel))
    if location[0]<0:
        gidx = "W"+str(xgid)
    else:
        gidx = "E"+str(xgid)
    if location[1]<0:
        gidy = "S"+ str(ygid)
    else:
        gidy="N"+str(ygid)
        
    GID=gidx+"-"+gidy
    return GID, xgid, ygid

def getrasteridx(x, y, w, h, x0, y0):
    idx = math.floor((x-x0)/w)
    idy = math.floor((y-y0)/h)
    return [int(idx),int(idy)]


def sample(lowerlocationsx, lowerlocationsy, number_of_samples, district_name, tiffile_1k, number_of_1_by_1=4, directory=''):
    
    tiffile_1k = "mean_output_1k.tif"
    crpdat = gdal.Open(tiffile_1k, gdal.GA_ReadOnly)
    crp_prob_np = np.asarray(crpdat.GetRasterBand(1).ReadAsArray())
    originX, pixelWidth, rx, originY, ry, pixelHeight = crpdat.GetGeoTransform()

    filenames = []
    for i in xrange(0, len(lowerlocationsx)):
        
        xycenter = (lowerlocationsx[i], lowerlocationsy[i])
        GID_center, xgid, ygid= GID_1k_func(xycenter) 
    
        res_pixel = 1000

        xycenter_GID = [xgid*res_pixel, ygid*res_pixel]
        xycenter_GID = [ -x if xycenter[m]<0 else x for m, x in enumerate(xycenter_GID)]

    ### specify grid resolution (grain, in m)
        grain = 100
        
        xoff =  xycenter_GID[0]- 5000 + grain/2
        yoff =  xycenter_GID[1]- 5000 + grain/2


        # xdim & ydim specify the number of cells to be sampled in x & y directions
        xdim = 100
        ydim = 100

        # specify grid resolution (grain, in m)
        grain = 100

        # generate the grid
        x, y = np.mgrid[yoff + ydim*(grain - 1):yoff - grain:-grain,
                        xoff:xoff + xdim*grain:grain]
        
        # write grid to the DataFrame with y and x columns
        grid = pd.DataFrame(np.vstack([x.ravel(), y.ravel()])).T
        grid.columns = ['y', 'x']

        # set up level ID's at the desired scales (res.pixel, in m)
        res_pixel = np.array([5000, 1000, grain])
        X2 = ((grid['x'] - xoff + 1.0/2*xdim)/res_pixel[0]).map(math.ceil)
        Y2 = ((grid['y'] - yoff + 1.0/2*ydim)/res_pixel[0]).map(math.ceil)
        L2 = cleanstrata(concatenate(X2, Y2, sep=''))

        X1 = ((X2*res_pixel[0] - (grid['x'] - xoff))/res_pixel[1]).map(math.ceil)
        Y1 = ((Y2*res_pixel[0] - (grid['y'] - yoff))/res_pixel[1]).map(math.ceil)
        L1 = cleanstrata(concatenate(X1, Y1, sep=''))

        X0 = ((X1*res_pixel[1] + (grid['x'] - xoff - X2*res_pixel[0]))/res_pixel[2]).map(math.ceil)
        Y0 = ((Y1*res_pixel[1] + (grid['y'] - yoff - Y2*res_pixel[0]))/res_pixel[2]).map(math.ceil)
        L0 = cleanstrata(concatenate(X0, Y0, sep=''))

        
        GID_1k =[GID_1k_func([y.ravel()[k], x.ravel()[k]])[0] for k in xrange(len(x.ravel()))]
        x_1k = [GID_1k_func([y.ravel()[k], x.ravel()[k]])[1]*1000 + 500 if  y.ravel()[k] >0 else -GID_1k_func([y.ravel()[k], x.ravel()[k]])[1]*1000 -500 for k in xrange(len(x.ravel()))]
        y_1k = [GID_1k_func([y.ravel()[k], x.ravel()[k]])[2]*1000 + 500 if  x.ravel()[k] >0 else -GID_1k_func([y.ravel()[k], x.ravel()[k]])[2]*1000 -500 for k in xrange(len(x.ravel()))]
      
        project = pyproj.Proj('+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs')
        
        xy_1k_latlon = [project(x_1k[k], y_1k[k], inverse=True) for k in xrange(len(x_1k))]
        xy_1k_latlon = zip(*xy_1k_latlon)
        
        GID_L0 = zip(GID_1k, L0)
        GID_100m = map(lambda (x, y): x+"-"+str(y-1), GID_L0)
        
        IDs = concatenate(L2, L1, L0, sep='')

        # update grid with level ID's
        grid['L2'] = L2
        grid['L1'] = L1
        grid['L0'] = map(lambda x: x-1, L0)
        grid['ID_unit'] = IDs
        grid['GID_1k'] = GID_1k
        grid['GID_100m'] = GID_100m
        grid['x_1k'] = list(xy_1k_latlon[0])
        grid['y_1k'] = list(xy_1k_latlon[1])
        
        # calculate selection probabilities for 1k-by-1k grids

        xy_1k = list(set(map(tuple,zip(*[x_1k, y_1k]))))
        xy_1k = zip(*map(lambda x: list(x), xy_1k))
        xy_1k.extend([tuple([pixelWidth]*len(x_1k)), tuple([pixelHeight]*len(x_1k)), tuple([originX]*len(x_1k)), tuple([originY]*len(x_1k))])
        crp_presence_idx = list(starmap(getrasteridx, zip(*xy_1k)))

        crp_prob_grid = [crp_prob_np[x[1]][x[0]] for x in crp_presence_idx]

        xy_1k = zip(*xy_1k[:2])
        xy_1k_pd = pd.DataFrame(xy_1k)
        xy_1k_pd.columns = ['x', 'y']
        X2 = ((xy_1k_pd['x'] - xoff + 1.0/2*xdim)/res_pixel[0]).map(math.ceil)
        Y2 = ((xy_1k_pd['y'] - yoff + 1.0/2*ydim)/res_pixel[0]).map(math.ceil)
        L2_1k = cleanstrata(concatenate(X2, Y2, sep=''))

        X1 = ((X2*res_pixel[0] - (xy_1k_pd['x'] - xoff))/res_pixel[1]).map(math.ceil)
        Y1 = ((Y2*res_pixel[0] - (xy_1k_pd['y'] - yoff))/res_pixel[1]).map(math.ceil)
        L1_1k = cleanstrata(concatenate(X1, Y1, sep=''))

        cutoff_1k = 0.3        
        select_prob_1k = map(lambda x: 0 if np.isnan(x) or x<cutoff_1k else x, crp_prob_grid)

        select_prob_1k =pd.DataFrame(zip(*[select_prob_1k, L2_1k, L1_1k]))
        select_prob_1k.columns =['prob', 'L2', 'L1']
        select_prob_1k['num_nonzero'] = select_prob_1k.prob>0
        
        size_1k_permit = select_prob_1k.groupby('L2').sum().num_nonzero.get_values()
        
        if any(size_1k_permit):
            samples = []
            L1_prob = []
            L2_prob = []

            size_1k = [int(x) if x< number_of_1_by_1 else number_of_1_by_1 for x in size_1k_permit]
            print(size_1k)
        # taking samples
            for (name, L2_group) in grid.groupby(['L2']):
            # take all L2 groups
                if(size_1k[name-1]>0):
                    grouped = L2_group.groupby(['L1'])
                    choices = [key for key in grouped.indices.keys()]
                    prob = select_prob_1k[select_prob_1k.L2==name].sort('L1').prob.get_values()
                    prob = map(lambda x: x*1.0/sum(prob), prob)
            # select L1 in each L2
                    sample = np.random.choice(choices, size_1k[name-1], p=prob)

                    for L1_value in sample:
                # in each selected L1 select 100x100 blocks
                        L1_group = grouped.get_group(L1_value)
                        choices = [j for j in xrange(len(L1_group.index))]
                        sample = random.sample(choices, number_of_samples)
                        for s in sample:
                            row = L1_group.irow(s)
                    # project back to longitude/latitude
                            row['x'], row['y'] = project(row['x'], row['y'], inverse=True)
                            samples.append(row)
                            L1_prob_tmp = prob[L1_value-1]
                            L1_prob.append(L1_prob_tmp)
                            L2_prob.append(number_of_samples/(1.0*res_pixel[1]/res_pixel[2])**2  * L1_prob_tmp)
                    
            df = pd.DataFrame(samples)
        # find probability of choosing L1
            df['prob_1k'] = L1_prob
        # find probability of choosing an 100x100 block
            df['prob'] = L2_prob

        
            if not os.path.exists('output/'+district_name):
                os.makedirs('output/'+district_name)

        
            if not os.path.exists('output/'+district_name+'/drone_flight_1k'):
                os.makedirs('output/'+district_name+'/drone_flight_1k')
        
            # write to csv
            filename = 'output/'+district_name+'/'+district_name+"_"+str(i)
            df.to_csv(filename + '.csv', index=False)

        # convert to kml
            spatial_csv_to_kml.csv_to_kml(filename + '.csv')

        # convert to gpx 
            cmd = 'gpsbabel -i kml -f "' +filename+'.kml" -o gpx -F "'+ filename+'_Waypoints.gpx"'
            c = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (stdout, stderr) = c.communicate()
            rc = c.returncode

            filenames.append(filename)
    return filenames



if __name__ == '__main__':
    # Example of running from command line:
    # ./sampling.py 38,3 10,8 40
    # the 1st argument is a list of longitudes, the 2nd is a list of latitudes and the 3rd is the number of samples
    print sample([float(x) for x in sys.argv[1].split(',')], [float(y) for y in sys.argv[2].split(',')], int(sys.argv[3]), str(sys.argv[4]), str(sys.argv[5]))
