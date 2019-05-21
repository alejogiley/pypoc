# -*- coding: utf-8 -*-

"""Main module."""

###########################################################################
#### imports                                                              #
###########################################################################

import sys, time

from modules import fmda
from modules import grdx

###########################################################################
#### functions/classes                                                    #
###########################################################################

def main(argv):
	
	# volume file
	volpath = "../examples/volume.dx"

	# execution time
	start = time.time()
	
	# run analysis
	data = fmda.mda2dict(argv)
		
	# format data dict
	grid, origin = grdx.dict2grd(data)
	
	# gaussian smooth
	kdde = grdx.grd4smooth(grid, 2.)
	
	# create dx object
	dtdx = grdx.grd2dx(kdde, origin)
	
	# save dx file
	dtdx.write(volpath) 
	
	# interpolate dx object
	grdx.dx4resample(volpath, 1) 
	
	# execution time
	end = time.time()
	
	# calculate time
	elp = end - start
	
	# print time
	print("%.2f s" % elp)

if __name__ == "__main__":
	main(sys.argv[1:])
