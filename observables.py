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
import mbii.lego_tools as util
import mbii.basic_simulation_info as info

#plt.style.use('y1a1')
def corrs_with_variance(self, data, niter=10, seed=9000, verbosity=1):
    """Compute the various 1h and 2h terms in the galaxy galaxy realspace correlation.
    Then compute bootstrap errorbars by subsampling."""
    bootstrap_iterations=[]

    res = groups._get_corrs(data.info, ctype=('s','c'), min_sep=44, max_sep=18e3, nbins=18, randoms=ran)

    np.random.seed(seed)
    if verbosity>0:
        print 'Constructing catalogue of random points'
    ran = util.construct_random_cat(data1.info,mask=None)

    for i in xrange(niter):
        if verbosity>0:
            print 'Iteration %d/%d'%(i+1,niter)
        indices = np.random.choice( int(len(data.info)), int(len(data.info)/(1.*niter)), replace=False )
        res_iter = groups._get_corrs(data.info[indices], ctype=('s','c'), min_sep=44, max_sep=18e3, nbins=18, randoms=ran)
        bootstrap_iterations.append(res_iter[1])

    return res[0], res[1], np.array(bootstrap_iterations)


def get_deltasigma(galaxy_data, particle_data,  min_sep=44, max_sep=2e3, binning='log', nbins=20, verbosity=1, downsampling_factor=1):

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

