import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
import copy
import mbii.lego_tools as util

def jackknife(data1, data2, options, verbosity=0):
	nsub = options['errors']['nsub']
	if verbosity>0:
		print 'Calculating jackknife errorbars - %dx%d subvolumes'%(nsub,nsub) 

	dx = 100./nsub
	if verbosity>0:
		print 'sub-box length : %3.3f h^-1 Mpc'%dx

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
					print data1['x'][mask1].mean()

				cat1 = treecorr.Catalog(x=data1['x'][mask1], y=data1['y'][mask1], z=data1['z'][mask1], a=data1['a1'][mask1], b=data1['a2'][mask1], c=data1['a3'][mask1])
				cat2 = treecorr.Catalog(x=data2['x'][mask2], y=data2['y'][mask2], z=data2['z'][mask2], a=data2['a1'][mask2], b=data2['a2'][mask2], c=data2['a3'][mask2])
				ed = treecorr.NVCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])

				ed.process(cat1,cat2)
				ED.append(copy.deepcopy(ed.xi))
				ed.clear()
				del(ed.xi)
				nprocessed+=1
				if verbosity>0:
					print '%d/%d'%(nprocessed,nsub*nsub*nsub)

	if verbosity>0:
		print 'Done subsampling.'
	return np.array(ED).std(axis=0)

