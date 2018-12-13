import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
import copy
#import mbii.lego_tools as util

period={'massiveblackii':100, 'illustris':75}

def jackknife(data1, data2, options, verbosity=0, nbins=6):
	nsub = options['errors']['nsub']
	if verbosity>0:
		print('Calculating jackknife errorbars - %dx%d subvolumes'%(nsub,nsub) )

	dx = period[options['simulation']]/nsub
	if verbosity>0:
		print('sub-box length : %3.3f h^-1 Mpc'%dx)

	vec=[]
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

				cat1 = treecorr.Catalog(x=data1['x'][mask1], y=data1['y'][mask1], z=data1['z'][mask1])
				cat2 = treecorr.Catalog(x=data2['x'][mask2], y=data2['y'][mask2], z=data2['z'][mask2])
				gg = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=nbins)

				gg.process(cat1,cat2)
				rcat11, rcat12, rcat21, rcat22 = randoms(cat1, cat2, period=period[options['simulation']])
				gg = finish_nn(gg, cat1, cat1, rcat11, rcat12, options, nbins)

				vec.append(copy.deepcopy(gg.xi))
				gg.clear()
				del(gg.xi)
				nprocessed+=1
				if verbosity>0:
					print('%d/%d'%(nprocessed,nsub*nsub*nsub))

	if verbosity>0:
		print('Done subsampling.')

	return np.array(vec).std(axis=0)

def randoms(cat1, cat2, period=100):
	# Initialise randoms
	# Fix the random seed so the randoms are always the same for a catalog of given length
	# Should probably migrate this to a config option

	# Need to fiddle with the range from which the randoms are drawn
	# to ensure they fall within the correct simulation volume
	np.random.seed(9000)

	rcat11 = treecorr.Catalog(x=np.random.rand(cat1.ntot)*period, y=np.random.rand(cat1.ntot)*period, z=np.random.rand(cat1.ntot)*period)
	rcat12 = treecorr.Catalog(x=np.random.rand(cat1.ntot)*period, y=np.random.rand(cat1.ntot)*period, z=np.random.rand(cat1.ntot)*period)
	rcat21 = treecorr.Catalog(x=np.random.rand(cat2.ntot)*period, y=np.random.rand(cat2.ntot)*period, z=np.random.rand(cat2.ntot)*period)
	rcat22 = treecorr.Catalog(x=np.random.rand(cat2.ntot)*period, y=np.random.rand(cat2.ntot)*period, z=np.random.rand(cat2.ntot)*period)

	return rcat11, rcat12, rcat21, rcat22

def finish_nn(corr, cat1, cat2, rcat1, rcat2, options, nbins):

	if cat2 is None:
		cat2 = cat1

	rr = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=nbins)
	rn = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=nbins)
	nr = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=nbins)

	# Do the pair counting
	print('Processing randoms',)
	print('RR',)
	rr.process(rcat1,rcat2)
	print('DR',)
	nr.process(cat1,rcat2)
	print('RD',)
	rn.process(rcat1,cat2)

	# Finish off
	xi, var = corr.calculateXi(rr, dr=nr, rd=rn)
	setattr(corr, 'xi', xi)

	return corr

