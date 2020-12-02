"""
* Sampling script that seeks to fit samples of certain sizes into boundaries and
*  then calculates mean and standard deviation of values from an underlying point dataset
*
* Used for Adeluyi et al. (2021)
*
* @author jonnyhuck
"""

from rtree import index
from pandas import DataFrame
from numpy.random import uniform
from statistics import mean, pstdev
from shapely.affinity import rotate
from shapely.geometry import Polygon
from geopandas import read_file, GeoSeries, GeoDataFrame

''' SETTINGS---------------------------------------------------------------- '''

sample_n = 100      # number of samples that you want
termination_f = 500 # number of iterations after which the algorithm gives up (*sample_n)
sample_sizes = [    # [min, max] size in hectares
    [0.5, 2],
    [2, 10],
    [0.01, 1],
    [1, 2],
    [2,5],
    [5,10]
    ]

''' ------------------------------------------------------------------------ '''

# update user
print("preparing data...")

# open dataset
data = read_file('./in/Yield_10m_point_seun.shp')

# get bounds
bounds = data.total_bounds

# initialise an rtree index
idx = index.Index()
for id, d in data.iterrows():
	idx.insert(id, d.geometry.bounds)

# read in the aoi
aoi = read_file("./in/olam_boundary3.shp")

# update user
print("calculating samples...")

# loop through each sample size
for sample_size in sample_sizes:

    # report current sample size
    print(sample_size)

    # lists to hold the results
    out_mean = []
    out_stdev = []
    out_geom = []

    # init/reset timer
    count = 0

    # keep generating samples until we have enough
    while len(out_geom) < sample_n:

        # give up if we reach the termination factor
        count +=1
        if count == sample_n * termination_f:
            print(f"{sample_size[0]}-{sample_size[1]} failed at {len(out_geom)} samples.")
            break

        # get the size of the random field
        size = (uniform(sample_size[0], sample_size[1]) * 10000)**0.5

        # random point to originate the box
        point = [int(uniform(bounds[0], bounds[2]-size)), int(uniform(bounds[1], bounds[3]-size))]

        # build the box from the point & rotate by 35 degrees to align with the dataset
        polygon = rotate(Polygon([(point[0], point[1]), (point[0], point[1]+size),
            (point[0]+size, point[1]+size), (point[0]+size, point[1])]), 45, origin='centroid')

        # test for inetrsection with the aoi (adapted for multiple polygons)
        intersects = False
        for i in range(len(aoi.index)):
            if polygon.within(aoi.geometry.iloc[i]):
                intersects = True

        # if no intersection, next sample
        if not intersects:
            continue

        # use spatial indexes to extract the possible then precise points in the sample field
        possible_matches = data.iloc[list(idx.intersection(polygon.bounds))]
        precise_matches = possible_matches[possible_matches.within(polygon)]

        # calculate stats and store in lists ready for output shapefile
        vals = list(precise_matches.GRID_CODE)
        out_mean.append(mean(vals))
        out_stdev.append(pstdev(vals))
        out_geom.append(polygon)

    # export result (if there is one to export)
    if len(out_geom) > 0:
        result = GeoDataFrame(DataFrame(data={'mean':out_mean, 'stdev':out_stdev}),
            geometry=GeoSeries(out_geom, crs=data.crs), crs=data.crs)
        result['area'] = result.geometry.area * 0.0001
        result.to_file(f"./out/{sample_size[0]}-{sample_size[1]}-samples.shp")

# report completed
print("done")
