# -*- Mode: python
#
# Licenses 
# Authors

"""\
Functions to analyze the simulation
===================================

Comments

"""

###########################################################################
#### imports                                                              #
###########################################################################

import numpy as np
import MDAnalysis as mda

from pypoc import TopologyGroup
from pypoc.clase import _Acyls, add_to_dict
from collections import defaultdict

###########################################################################
#### functions / classses                                                 #
###########################################################################

# !!!!!!!!! Paralellize
def order2rank(u):
	"""Comments"""
	
	# dict of resid, 
	# bond, and angles
	data = {}

	# load trajectory
	u = mda.Universe(argv[0], argv[1]) 
	# lipids residues
	lipids = u.select_atoms("resname POPE POPG")
	
	# acyl atoms
	"""---- POPC DLPC ALL ATOMS
	name  = ["C2%i"%i for j in (range(1,10), range(10,19)) for i in j]
	name += ["C3%i"%i for j in (range(1,10), range(10,17)) for i in j]
	"""
	name  = ["C%iB" %i for i in (range(1,5))]
	name += ["C%iA" %i for i in (1,3,4)]
	name += ["D2A"]
	lips  = lipids.select_atoms("name " + " ".join(name), updating=True)
	
	# topology dict
	ad  = lips.bonds.topDict
	tg  = ad[('C1', 'C1')]
	tg += ad[('C1', 'C3')]
	
	"""---- ALL ATOMS
	tg  = ad[('CL', 'CTL2')]
	tg += ad[('CTL2', 'CTL2')]
	tg += ad[('CTL2', 'CTL3')]
	tg += ad[('CTL2', 'CEL1')]
	tg += ad[('CEL1', 'CEL1')]
	"""

	for ts in u.trajectory:
		############################################################
		### integer coordinates for each atom1			 ###
		############################################################
		# get coordinates
		x = list(zip(*tg.atom1.positions))[0] 
		y = list(zip(*tg.atom1.positions))[1]
		z = list(zip(*tg.atom1.positions))[2]
		# round positions
		x = map(round, x)
		y = map(round, y)
		z = map(round, z)
		# data dict keys
		ids = list(map(lambda a,b,c: 
			",".join(map(str,[a,b,c])), x, y, z))
                ############################################################
                ### dict of bond angles for each residue		 ###
                ############################################################	
		# get resid & angles
		res = tg.atom1.resids
		ang = tg._srankSlow(pbc=False)
		# save in dict
		resang = defaultdict(list)
		for key, val in zip(res,ang):
			resang[key].append(val)
                ############################################################
                ### calculate z-angles for each residue 	         ###
                ############################################################
		# length of dict is one per residue
		angles = [_Acyls.average(resang[k]) 
				for k in list(set(res))]
		# length of list one per bond
		series = [angles[i] for i in res-res[0]]
		############################################################
                ### save values in data dict 	                         ###
                ############################################################
		add_to_dict(data, ids, tg.atom1.resids, tg.atom1.ix, series)
	
	return data

