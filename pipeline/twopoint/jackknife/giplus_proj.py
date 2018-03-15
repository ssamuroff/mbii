import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
import copy
import mbii.lego_tools as util
from halotools.mock_observables.alignments import gi_plus_projected

def jackknife(data1, data2, options, verbosity=0, rpbins=None):
	nsub = options['errors']['nsub']
	if verbosity>0:
		print 'Calculating jackknife errorbars - %dx%d subvolumes'%(nsub,nsub) 

	dx = 100./nsub
	if verbosity>0:
		print 'sub-box length : %3.3f h^-1 Mpc'%dx

	GI=[]
	nprocessed=0

	for i in xrange(nsub-1):
		# x axis box
		xmask1 = (data1['x']<dx*i) | (data1['x']>dx*(i+1))
		xmask2 = (data2['x']<dx*i) | (data2['x']>dx*(i+1))

		for j in xrange(nsub):
			# y axis box
			ymask1 = (data1['y']<dx*j) | (data1['y']>dx*(j+1))
			ymask2 = (data2['y']<dx*j) | (data2['y']>dx*(j+1))

			for k in xrange(nsub):
				# z axis box
				zmask1 = (data1['z']<dx*k) | (data1['z']>dx*(k+1))
				zmask2 = (data2['z']<dx*k) | (data2['z']>dx*(k+1))

				# Combined mask for 3D subvolume
				mask1 = xmask1 & ymask1 & zmask1
				mask2 = xmask2 & ymask2 & zmask2

				if verbosity>0:
					print data1['x'][mask1].mean()

				cat1 = data1[mask1]
				cat2 = data2[mask2]
				gi = compute_giplus(cat1, cat2, options, period=100., rpbins=rpbins)

				GI.append(copy.deepcopy(gi))
				nprocessed+=1
				if verbosity>0:
					print '%d/%d'%(nprocessed,nsub*nsub*nsub)

	if verbosity>0:
		print 'Done subsampling.'
	return np.array(GI).std(axis=0)

def compute_giplus(cat1, cat2, options, period=100., rpbins=None):
	avec = np.vstack((cat2['a1'], cat2['a2'])).T
	pvec1 = np.vstack((cat1['x'], cat1['y'], cat1['z'])).T
	pvec2 = np.vstack((cat2['x'], cat2['y'], cat2['z'])).T
	evec = np.sqrt(cat2['e1']*cat2['e1'] + cat2['e2']*cat2['e2'])

	if (options['2pt']['binning']=='log') and (rpbins is None):
		rpbins = np.logspace(np.log10(options['2pt']['rpmin']), np.log10(options['2pt']['rpmax']), options['2pt']['nrpbin'])
	elif (options['2pt']['binning']=='equal') and (rpbins is None):
		rpbins = util.equalise_binning(options['2pt']['rpmin'], options['2pt']['rpmax'], options['2pt']['nbin'])
	pi_max = options['2pt']['pi_max']

	gip = gi_plus_projected(pvec2, avec, evec, pvec1, rpbins, pi_max, period=period, num_threads=1) 

	return gip
