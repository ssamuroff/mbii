import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
import copy
#import mbii.lego_tools as util
from halotools.mock_observables.alignments import gi_plus_projected

period={'massiveblackii':100, 'illustris':75}

def jackknife(data1, data2, options, verbosity=0, rpbins=None, nbins=6):
	nsub = options['errors']['nsub']
	if verbosity>0:
		print('Calculating jackknife errorbars - %dx%d subvolumes'%(nsub,nsub) )

	dx = period[options['simulation']]/nsub
	if verbosity>0:
		print('sub-box length : %3.3f h^-1 Mpc'%dx)

	GI=[]
	nprocessed=0

	for i in range(nsub-1):
		# x axis box
		xmask1 = (data1['x']>dx*i) & (data1['x']<dx*(i+1))
		xmask2 = (data2['x']>dx*i) & (data2['x']<dx*(i+1))

		for j in range(nsub):
			# y axis box
			ymask1 = (data1['y']>dx*j) & (data1['y']<dx*(j+1))
			ymask2 = (data2['y']>dx*j) & (data2['y']<dx*(j+1))

			for k in range(nsub):
				# z axis box
				zmask1 = (data1['z']>dx*k) & (data1['z']<dx*(k+1))
				zmask2 = (data2['z']>dx*k) & (data2['z']<dx*(k+1))

				# Combined mask for 3D subvolume
				mask1 = np.invert(xmask1 & ymask1 & zmask1)
				mask2 = np.invert(xmask2 & ymask2 & zmask2)

				if verbosity>0:
					print(data1['x'][mask1].mean())

				cat1 = data1[mask1]
				cat2 = data2[mask2]
				gi = compute_giplus(cat1, cat2, options, rbins=rpbins, period=period[options['simulation']], nbins=nbins)

				GI.append(copy.deepcopy(gi))
				nprocessed+=1
				if verbosity>0:
					print('%d/%d'%(nprocessed,nsub*nsub*nsub))

	if verbosity>0:
		print('Done subsampling.')
	gi0 = np.mean(GI, axis=0)
	R2 = np.sum([(f - gi0)*(f - gi0) for f in GI], axis=0)
	coeff = (nsub**3 - 1.)/nsub**3
	dGI = np.sqrt(coeff * R2)

	return dGI

def compute_giplus(cat1, cat2, options, period=100., rbins=None, nbins=6):

	aname = 'a%d'
	ename = 'e%d'
	if ('shapes_suffix' in options['2pt'].keys()):
		aname+=options['2pt']['shapes_suffix']
		ename+=options['2pt']['shapes_suffix']

	avec = np.vstack((cat2[aname%1], cat2[aname%2])).T
	pvec1 = np.vstack((cat1['x'], cat1['y'], cat1['z'])).T
	pvec2 = np.vstack((cat2['x'], cat2['y'], cat2['z'])).T
	evec = np.sqrt(cat2[ename%1]*cat2[ename%1] + cat2[ename%2]*cat2[ename%2])

	if (options['2pt']['binning']=='log') and (rbins is None):
		rbins = np.logspace(np.log10(options['2pt']['rmin']), np.log10(options['2pt']['rmax']), nbins+1)
	elif (options['2pt']['binning']=='equal') and (rbins is None):
		rbins = util.equalise_binning(options['2pt']['rmin'], options['2pt']['rmax'], nbins+1)
	pi_max = options['2pt']['pi_max']

	mask2=avec.T[0]!=0.0

	gip = gi_plus_projected(pvec2[mask2], avec[mask2], evec[mask2], pvec1, rbins, pi_max, period=period, num_threads=1) 

	return gip
