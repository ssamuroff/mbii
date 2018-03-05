import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
import copy
from halotools.mock_observables.alignments import gi_plus_3d 

def jackknife(data1, data2, options, verbosity=0):
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
				gi = compute_giplus(cat1, cat2, options, period=100.)

				GI.append(copy.deepcopy(gi))
				nprocessed+=1
				if verbosity>0:
					print '%d/%d'%(nprocessed,nsub*nsub*nsub)

	if verbosity>0:
		print 'Done subsampling.'
	return np.array(GI).std(axis=0)

def compute_giplus(cat1, cat2, options, period=100.):
	avec = np.vstack((cat2['a1'], cat2['a2'], cat2['a3'])).T
	pvec1 = np.vstack((cat1['x'], cat1['y'], cat1['z'])).T
	pvec2 = np.vstack((cat2['x'], cat2['y'], cat2['z'])).T
	evec = np.sqrt(cat2['e1']*cat2['e1'] + cat2['e2']*cat2['e2'])

	rbins = np.logspace(np.log10(options['2pt']['rmin']), np.log10(options['2pt']['rmax']), options['2pt']['nbin'] )

	gip = gi_plus_3d(pvec2, avec, evec, pvec1, rbins, period=period, num_threads=1) 

	return gip
