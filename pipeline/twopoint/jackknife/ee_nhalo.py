import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
import copy
from halotools.mock_observables.alignments import ee_3d_one_two_halo_decomp as ee_3d
#import mbii.lego_tools as util


def jackknife(data1, data2, options, verbosity=0, rbins=None):
	nsub = options['errors']['nsub']
	if verbosity>0:
		print 'Calculating jackknife errorbars - %dx%d subvolumes'%(nsub,nsub) 

	dx = 100./nsub
	if verbosity>0:
		print 'sub-box length : %3.3f h^-1 Mpc'%dx

	EE1=[]
	EE2=[]
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
				ee_1h, ee_2h = compute_ee(cat1, cat2, options, period=100., rbins=rbins)

				EE1.append(copy.deepcopy(ee_1h))
				EE2.append(copy.deepcopy(ee_2h))
				nprocessed+=1
				if verbosity>0:
					print '%d/%d'%(nprocessed,nsub*nsub*nsub)

	if verbosity>0:
		print 'Done subsampling.'

	return np.array(EE1).std(axis=0), np.array(EE2).std(axis=0)

def compute_ee(cat1, cat2, options, period=100., rbins=None):
	avec1 = np.vstack((cat1['a1'], cat1['a2'], cat1['a3'])).T
	avec2 = np.vstack((cat2['a1'], cat2['a2'], cat2['a3'])).T
	pvec1 = np.vstack((cat1['x'], cat1['y'], cat1['z'])).T
	pvec2 = np.vstack((cat2['x'], cat2['y'], cat2['z'])).T
	evec = np.sqrt(cat2['e1']*cat2['e1'] + cat2['e2']*cat2['e2'])

	if (options['2pt']['binning']=='log') and (rbins is None):
		rbins = np.logspace(np.log10(options['2pt']['rmin']), np.log10(options['2pt']['rmax']), options['2pt']['nbin'])
	elif (options['2pt']['binning']=='equal') and (rbins is None):
		rbins = util.equalise_binning(options['2pt']['rmin'], options['2pt']['rmax'], options['2pt']['nbin'])

	hids1 = cat1['halo_id']
	hids2 = cat2['halo_id']

	ee_1h, ee_2h = ee_3d(pvec2, avec2, hids2, pvec1, avec1, hids1, rbins, period=period, num_threads=1) 

	return ee_1h, ee_2h
