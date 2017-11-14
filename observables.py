#import MySQLdb as mdb
import pymysql as mdb
import numpy as np
from numpy.core.records import fromarrays
import pylab as plt
import treecorr
import tools.plots as pl
plt.switch_backend('pdf')
import halotools as ht
import halotools.mock_observables as pretending
import mbii.lego_tools as tools
import mbii.basic_simulation_info as info

#plt.style.use('y1a1')



def get_deltasigma(self, galaxy_data, particle_data,  min_sep=44, max_sep=2e3, binning='log', nbins=20, verbosity=1, downsampling_factor=1):

	    if verbosity>0:
    	print 'Will construct %s - %s \Delta \Sigma profile'%ctype

    	# Decide on an appropriate binning scheme
    	if (binning.lower()=='log'):
    		rbins = np.logspace(np.log10(min_sep), np.log10(max_sep), nbins )
    	elif (binning.lower()=='linear'):
    		rbins = np.linspace(min_sep, max_sep, nbins )

    	if verbosity>1:
    		print 'Will use %s binning:'%binning, rbins

    	# Parse the mask
    	#mask1 = tools.choose_cs_mask(data,ctype[0])
    	#mask2 = tools.choose_cs_mask(data,ctype[1])

    	pos1 = pretending.return_xyz_formatted_array(particle_data['x'], particle_data['y'], particle_data['z']) #, mask = mask1)
    	pos2 = pretending.return_xyz_formatted_array(galaxy_data['x'], galaxy_data['y'], galaxy_data['z']) #, mask = mask2)

    	R = np.sqrt(np.array(rbins)[1:]*np.array(rbins)[:-1]) 

    	rp, DeltaSigma = pretending.delta_sigma(pos2, pos1, particle_data['mass'], downsampling_factor,
                    rbins, info.Lbox, cosmology=info.cosmology, num_threads='max')

    	return rp, DeltaSigma

