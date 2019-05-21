"""Example"""

import time
import popypc
import MDAnalysis as mda

# execution time
start = time.time()

# ----------------------

# read files
u = mda.Universe('data/clc-e1_pc.tpr',
		 'data/clc-e1_pc.pdb')

# ----------------------

# run analysis
data = analysis(u)

# ----------------------

# format data dict
grid, origin = dat2grd(data)
# gaussian smooth
kdde = smooth(grid, 1.0)
# create dx object
dtdx = grd2dx(kdde, origin)
# save dx file
dtdx.write("volume.dx")	
# interpolate dx object
dxgrid("volume.dx", 6) 

# ----------------------

# execution time
end = time.time()
# calculate time
elp = end - start
# print time
print("%.2f s" % elp)
