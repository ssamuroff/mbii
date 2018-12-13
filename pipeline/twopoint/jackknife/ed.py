import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
import copy
from halotools.mock_observables.alignments import ed_3d
#import mbii.lego_tools as util

period={'massiveblackii':100, 'illustris':75}

def jackknife(data1, data2, options, verbosity=0, rbins=None, nbins=6):
	nsub = options['errors']['nsub']
	if verbosity>0:
		print('Calculating jackknife errorbars - %dx%d subvolumes'%(nsub,nsub) )

	dx = period[options['simulation']]/nsub
	if verbosity>0:
		print('sub-box length : %3.3f h^-1 Mpc'%dx)

	ED=[]
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
				ed = compute_ed(cat1, cat2, options, period=period[options['simulation']], rbins=rbins, nbins=nbins)

				ED.append(copy.deepcopy(ed))
				nprocessed+=1
				if verbosity>0:
					print('%d/%d'%(nprocessed,nsub*nsub*nsub))

	if verbosity>0:
		print('Done subsampling.')

	ed0 = np.mean(ED, axis=0)
	R2 = np.sum([(f - ed0)*(f - ed0) for f in ED], axis=0)
	coeff = (nsub**3 - 1.)/nsub**3
	dED = np.sqrt(coeff * R2)

	return dED

	#np.array(ED).std(axis=0)

def compute_ed(cat1, cat2, options, period=100., rbins=None, nbins=6):

	aname = 'a%d'
	ename = 'e%d'
	if ('shapes_suffix' in options['2pt'].keys()):
		aname+=options['2pt']['shapes_suffix']
		ename+=options['2pt']['shapes_suffix']

	avec1 = np.vstack((cat1[aname%1], cat1[aname%2], cat1[aname%3])).T
	avec2 = np.vstack((cat2[aname%1], cat2[aname%2], cat2[aname%3])).T
	pvec1 = np.vstack((cat1['x'], cat1['y'], cat1['z'])).T
	pvec2 = np.vstack((cat2['x'], cat2['y'], cat2['z'])).T
	evec = np.sqrt(cat2[ename%1]*cat2[ename%1] + cat2[ename%2]*cat2[ename%2])

	if (options['2pt']['binning']=='log') and (rbins is None):
		rbins = np.logspace(np.log10(options['2pt']['rmin']), np.log10(options['2pt']['rmax']), nbins+1)
	elif (options['2pt']['binning']=='equal') and (rbins is None):
		rbins = util.equalise_binning(options['2pt']['rmin'], options['2pt']['rmax'], nbins+1)

	mask1=avec1.T[0]!=0.0
	mask2=avec2.T[0]!=0.0

	ed = ed_3d(pvec2[mask2], avec2[mask2], pvec1[mask1], rbins, period=period, num_threads=1) 

	return ed

def jackknife_treecorr(data1, data2, options, verbosity=0):
	nsub = options['errors']['nsub']
	if verbosity>0:
		print('Calculating jackknife errorbars - %dx%d subvolumes'%(nsub,nsub) )

	dx = 100./nsub
	if verbosity>0:
		print('sub-box length : %3.3f h^-1 Mpc'%dx)

	ED=[]
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
					print(data1['x'][mask1].mean())

				cat1 = treecorr.Catalog(x=data1['x'][mask1], y=data1['y'][mask1], z=data1['z'][mask1], a=data1['a1'][mask1], b=data1['a2'][mask1], c=data1['a3'][mask1])
				cat2 = treecorr.Catalog(x=data2['x'][mask2], y=data2['y'][mask2], z=data2['z'][mask2], a=data2['a1'][mask2], b=data2['a2'][mask2], c=data2['a3'][mask2])
				ed = treecorr.NVCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])

				ed.process(cat1,cat2)
				ED.append(copy.deepcopy(ed.xi))
				ed.clear()
				del(ed.xi)
				nprocessed+=1
				if verbosity>0:
					print('%d/%d'%(nprocessed,nsub*nsub*nsub))

	if verbosity>0:
		print('Done subsampling.')
	return np.array(ED).std(axis=0)

