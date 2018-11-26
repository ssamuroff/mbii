import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
import copy
from halotools.mock_observables.alignments import ed_3d
from halotools.mock_observables.alignments import ee_3d
#from halotools.mock_observables.two_point_clustering.tpcf import tpcf as gg_3d 
from halotools.mock_observables.alignments import gi_plus_projected
from halotools.mock_observables.alignments import ii_plus_projected
#import mbii.lego_tools as util
import scipy.optimize as opt


def func(x, a, b):
	return 10**b * x**a


def gg_3d(pvec2, pvec1, rbins, rvec1, rvec2, options, nbins=6):
	# Turn the arrays into TreeCorr readable catalogues
	cat1 = treecorr.Catalog(x=pvec1.T[0], y=pvec1.T[1], z=pvec1.T[2])
	cat2 = treecorr.Catalog(x=pvec2.T[0], y=pvec2.T[1], z=pvec2.T[2])

	# Initialise the correlation function
	gg = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=nbins)
	# Do the calculation
	gg.process(cat1,cat2)

	# Process the randoms
	rcat1 = treecorr.Catalog(x=rvec1.T[0], y=rvec1.T[1], z=rvec1.T[2])
	rcat2 = treecorr.Catalog(x=rvec2.T[0], y=rvec2.T[1], z=rvec2.T[2])
	gg = finish_nn(gg, cat1, cat2, rcat1, rcat2, options, nbins=nbins)

	R = np.exp(gg.logr)
	xi = gg.xi

	return xi

def finish_nn(corr, cat1, cat2, rcat1, rcat2, options, nbins=6):

	if cat2 is None:
		cat2 = cat1

	rr = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=nbins)
	rn = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=nbins)
	nr = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=nbins)

	# Do the pair counting
	#print('Processing randoms',)
	#print('RR',)
	rr.process(rcat1,rcat2)
	#print('DR',)
	nr.process(cat1,rcat2)
	#print('RD',)
	rn.process(rcat1,cat2)

	# Finish off
	xi, var = corr.calculateXi(rr, dr=nr, rd=rn)
	setattr(corr, 'xi', xi)

	return corr


def get_randoms(npts, period, bounds=None):
	out = np.zeros(npts, dtype=[('x', float), ('y', float), ('z', float)])
	out['x'] = np.random.rand(npts) * period
	out['y'] = np.random.rand(npts) * period
	out['z'] = np.random.rand(npts) * period

	if bounds is not None:
		for i in range(8000):
			mask = (out['x']>bounds[0][0]) & (out['x']<bounds[0][1]) & (out['y']>bounds[1][0]) & (out['y']<bounds[1][1]) & (out['z']>bounds[2][0]) & (out['z']<bounds[2][1])
			ngal = out['x'][mask].size
			if ngal==0:
				break
			out['x'][mask] = np.random.rand(ngal) * period
			out['y'][mask] = np.random.rand(ngal) * period
			out['z'][mask] = np.random.rand(ngal) * period

	return out

period={'massiveblackii':100, 'illustris':75}
measurement_functions = {'gg': gg_3d, 'ed': ed_3d, 'ee': ee_3d, 'gi_plus_projected': gi_plus_projected, 'ii_plus_projected': ii_plus_projected}

def jackknife(correlations, data1, data2, options, verbosity=0, nbins=[6]*5):
	nsub = options['errors']['nsub']

	if verbosity>0:
		print ('Calculating jackknife errorbars - %dx%d subvolumes'%(nsub,nsub) )

	#fn = measurement_functions[correlation]

	dx = period[options['simulation']]/nsub
	if verbosity>0:
		print( 'sub-box length : %3.3f h^-1 Mpc'%dx)

	F=[]
	nprocessed=0

	randoms1 = get_randoms(data1.size, period[options['simulation']])
	randoms2 = get_randoms(data2.size, period[options['simulation']])

	for i in range(nsub):
		# x axis box
		xmask1 = (data1['x']<dx*i) | (data1['x']>dx*(i+1))
		xmask2 = (data2['x']<dx*i) | (data2['x']>dx*(i+1))
		xmask1_randoms = (randoms1['x']<dx*i) | (randoms1['x']>dx*(i+1))
		xmask2_randoms = (randoms2['x']<dx*i) | (randoms2['x']>dx*(i+1))


		for j in range(nsub):
			# y axis box
			ymask1 = (data1['y']<dx*j) | (data1['y']>dx*(j+1))
			ymask2 = (data2['y']<dx*j) | (data2['y']>dx*(j+1))
			ymask1_randoms = (randoms1['y']<dx*j) | (randoms1['y']>dx*(j+1))
			ymask2_randoms = (randoms2['y']<dx*j) | (randoms2['y']>dx*(j+1))


			for k in range(nsub):
				# z axis box
				zmask1 = (data1['z']<dx*k) | (data1['z']>dx*(k+1))
				zmask2 = (data2['z']<dx*k) | (data2['z']>dx*(k+1))
				zmask1_randoms = (randoms1['z']<dx*k) | (randoms1['z']>dx*(k+1))
				zmask2_randoms = (randoms2['z']<dx*k) | (randoms2['z']>dx*(k+1))


				# Combined mask for 3D subvolume
				mask1 = xmask1 & ymask1 & zmask1
				mask2 = xmask2 & ymask2 & zmask2
				rmask1 = xmask1_randoms & ymask1_randoms & zmask1_randoms
				rmask2 = xmask2_randoms & ymask2_randoms & zmask2_randoms


				cat1 = data1[mask1]
				cat2 = data2[mask2]

				dd = [compute(c, cat1, cat2, options, randoms1[rmask1], randoms2[rmask2], nbins=n) for a,n in zip(correlations,nbins) ]
				F.append(np.concatenate(dd))

				nprocessed+=1
				if verbosity>0:
					print ('%d/%d'%(nprocessed,nsub*nsub*nsub))

	if verbosity>0:
		print( 'Done subsampling.')

	f0 = np.mean(F, axis=0)
	cov = np.zeros((len(f0),len(f0)))

	for f1,m1 in zip(F,f0):
		for f2,m2 in zip(F,f0):
			import pdb ; pdb.set_trace()





	R2 = np.sum([(f - f0)*(f - f0) for f in F], axis=0)
	coeff = (nsub**3 - 1.)/nsub**3
	cov = np.sqrt(coeff * cov)

	#import pdb ; pdb.set_trace()

	return dF, dR

def compute(correlation, cat1, cat2, options, randoms1, randoms2, rbins=None, nbins=6):
	avec1 = np.vstack((cat1['a1'], cat1['a2'], cat1['a3'])).T
	avec2 = np.vstack((cat2['a1'], cat2['a2'], cat2['a3'])).T
	avec1_2d = np.vstack((cat1['a1'], cat1['a2'])).T
	avec2_2d = np.vstack((cat2['a1'], cat2['a2'])).T
	pvec1 = np.vstack((cat1['x'], cat1['y'], cat1['z'])).T
	pvec2 = np.vstack((cat2['x'], cat2['y'], cat2['z'])).T
	evec1 = np.sqrt(cat1['e1']*cat1['e1'] + cat1['e2']*cat1['e2'])
	evec2 = np.sqrt(cat2['e1']*cat2['e1'] + cat2['e2']*cat2['e2'])

	rvec1 = np.vstack((randoms1['x'], randoms1['y'], randoms1['z'])).T
	rvec2 = np.vstack((randoms2['x'], randoms2['y'], randoms2['z'])).T

	if (options['2pt']['binning']=='log') and (rbins is None):
		rbins = np.logspace(np.log10(options['2pt']['rmin']), np.log10(options['2pt']['rmax']), nbins+1)
	elif (options['2pt']['binning']=='equal') and (rbins is None):
		rbins = util.equalise_binning(options['2pt']['rmin'], options['2pt']['rmax'], nbins+1)
		

	pi_max = options['2pt']['pi_max']
	mask1 = (avec1.T[0]!=0.0) & (pvec1.T[0]!=0.0) 
	mask2 = (avec2.T[0]!=0.0) & (pvec2.T[0]!=0.0) 

	fn = measurement_functions[correlation]
	if (correlation.lower()=='ed'):
		return fn(pvec2[mask2], avec2[mask2], pvec1[mask1], rbins, num_threads=1) 
	elif (correlation.lower()=='ee'):
		return fn(pvec2[mask2], avec2[mask2], pvec1[mask1], avec1[mask1], rbins, num_threads=1) 
	elif (correlation.lower()=='gi_plus_projected'):
		return fn(pvec2[mask2], avec2_2d[mask2], evec2[mask2], pvec1, rbins, pi_max, randoms1=rvec1, randoms2=rvec2, num_threads=1) 
	elif (correlation.lower()=='ii_plus_projected'):
		return fn(pvec2[mask2], avec2_2d[mask2], evec2[mask2], pvec1[mask1], avec1_2d[mask1], evec1[mask1], rbins, pi_max, randoms1=rvec1, randoms2=rvec2, num_threads=1) 
	elif (correlation.lower()=='gg'):
		return fn(pvec2[mask2], pvec1[mask1], rbins, rvec2, rvec1, options, nbins=nbins) 
	else:
		raise ValueError('Unknown correlation function: %s.'%correlation)

def jackknife_treecorr(data1, data2, options, verbosity=0):
	nsub = options['errors']['nsub']
	if verbosity>0:
		print( 'Calculating jackknife errorbars - %dx%d subvolumes'%(nsub,nsub) )

	dx = 100./nsub
	if verbosity>0:
		print( 'sub-box length : %3.3f h^-1 Mpc'%dx)

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
					print( data1['x'][mask1].mean())

				cat1 = treecorr.Catalog(x=data1['x'][mask1], y=data1['y'][mask1], z=data1['z'][mask1], a=data1['a1'][mask1], b=data1['a2'][mask1], c=data1['a3'][mask1])
				cat2 = treecorr.Catalog(x=data2['x'][mask2], y=data2['y'][mask2], z=data2['z'][mask2], a=data2['a1'][mask2], b=data2['a2'][mask2], c=data2['a3'][mask2])
				ed = treecorr.NVCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])

				ed.process(cat1,cat2)
				ED.append(copy.deepcopy(ed.xi))
				ed.clear()
				del(ed.xi)
				nprocessed+=1
				if verbosity>0:
					print( '%d/%d'%(nprocessed,nsub*nsub*nsub))

	if verbosity>0:
		print( 'Done subsampling.')
	return np.array(ED).std(axis=0)

