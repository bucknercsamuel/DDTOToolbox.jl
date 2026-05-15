from osgeo import gdal
import cv2
import os
import numpy as np
import pickle

# # ..:: MSL Easy Map ::..
# # Extract coordinates from UE4 manually
# path_mapfolder = os.getcwd() + '\\AirSim\\maps\\msl_test_easy\\'
# path_heightmap = path_mapfolder + 'height_map.png'
# minX_UE4 = -84000
# maxX_UE4 = +105000
# minY_UE4 = -69000
# maxY_UE4 = +107000
# minZ_UE4 = +5200
# maxZ_UE4 = +13230
# nedX_UE4 = +6370
# nedY_UE4 = +16660
# nedZ_UE4 = +8570

# # ..:: MSL Hard Map ::..
# # Extract coordinates from UE4 manually
# path_mapfolder = os.getcwd() + '\\AirSim\\maps\\msl_test_hard\\'
# path_heightmap = path_mapfolder + 'height_map.png'
# minX_UE4 = -73800
# maxX_UE4 = +77400
# minY_UE4 = -81000
# maxY_UE4 = +70200
# minZ_UE4 = -16400
# maxZ_UE4 = +33050
# nedX_UE4 = +1608
# nedY_UE4 = -31840
# nedZ_UE4 = -6518

# ..:: Plateau Medium Map ::..
# Extract coordinates from UE4 manually
path_mapfolder = os.getcwd() + '\\AirSim\\maps\\plateau_test_medium\\'
path_heightmap = path_mapfolder + 'height_map.png'
minX_UE4 = -61400
maxX_UE4 = +77200
minY_UE4 = -72400
maxY_UE4 = +66200
minZ_UE4 = -5248
maxZ_UE4 = +3672
nedX_UE4 = +9900
nedY_UE4 = -16140
nedZ_UE4 = -1308

# # ..:: Dunes Hard Map ::..
# # Extract coordinates from UE4 manually
# path_mapfolder = os.getcwd() + '\\AirSim\\maps\\dunes_test_hard\\'
# path_heightmap = path_mapfolder + 'height_map.png'
# minX_UE4 = -86400
# maxX_UE4 = +77400
# minY_UE4 = -84000
# maxY_UE4 = +79800
# minZ_UE4 = -19158
# maxZ_UE4 = +11972
# nedX_UE4 = +3470
# nedY_UE4 = +34710
# nedZ_UE4 = -9728


# ..:: Processing ::..
nedX_map_fraction = (nedX_UE4 - minX_UE4)/(maxX_UE4 - minX_UE4)
nedY_map_fraction = (nedY_UE4 - minY_UE4)/(maxY_UE4 - minY_UE4)
nedZ_map_fraction = (nedZ_UE4 - minZ_UE4)/(maxZ_UE4 - minZ_UE4)

# Obtain height map and data
raster = gdal.Open(path_heightmap)
array = raster.ReadAsArray()
min_array = array.min()
max_array = array.max()
norm_array = array - min_array
norm_array = norm_array/(norm_array.max())
norm_array = norm_array*65535
norm_array = norm_array.astype('uint16')
map_widthX = array.shape[1]
map_widthY = array.shape[0]
map_widthZ = 0.01*(maxZ_UE4 - minZ_UE4) # Convert cm to m

# Find min/max NED coordinates across map
minX_ned = -map_widthX*nedX_map_fraction
minY_ned = -map_widthY*nedY_map_fraction
minZ_ned = -map_widthZ*nedZ_map_fraction
maxX_ned = minX_ned + map_widthX
maxY_ned = minY_ned + map_widthY
maxZ_ned = minZ_ned + map_widthZ

# Create lookup table
lookup_table = {}
for j,x in enumerate(range(int(minX_ned), int(maxX_ned))):
    for i,y in enumerate(range(int(minY_ned), int(maxY_ned))):
        lookup_table[(x,y)] = norm_array[i][j]*(map_widthZ/65535.) + minZ_ned

# Save lookup table
print("NED map extrema coordinates: ")
print("Min-X NED: ", minX_ned)
print("Max-X NED: ", maxX_ned)
print("Min-Y NED: ", minY_ned)
print("Max-Y NED: ", maxY_ned)
f = open(path_mapfolder + 'lookup_table.pkl', 'wb')
pickle.dump(lookup_table, f)
f.close()