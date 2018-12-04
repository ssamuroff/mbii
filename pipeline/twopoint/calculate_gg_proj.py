import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
#import mbii.lego_tools as util
from mbii.pipeline.twopoint.jackknife import gg_proj as errors 
from halotools.mock_observables.two_point_clustering import wp  

periods={'massiveblackii':100, 'illustris':75}

def compute0(options, binning):
	
	if 'positions' in options['2pt'].keys():
		posfile = options['2pt']['positions']
		data = fi.FITS(posfile)[-1].read()
	else:
		posfile = options['2pt']['shapes']
		data = fi.FITS(posfile)[-1].read()

	print('Position data : %s'%posfile)

	splitflag=options['2pt']['split']
	if 'npart_cut' in options['2pt'].keys():
		nlow = options['2pt']['npart_cut']
		print('Imposing additional cut at npart>%d'%nlow)
		data = data[(data['npart_dm']>nlow)]
	else:
		nlow=-1

	if 'npart_cut_upper' in options['2pt'].keys():
		nhigh = options['2pt']['npart_cut_upper']
		print('Imposing additional cut at npart<%d'%nhigh)
		data = data[(data['npart_dm']<nhigh)]
	else:
		nhigh = -1

	if splitflag is not None:
		name = options['2pt']['split']
		print('Dividing catalogue by %s'%name)
		mask = (data[name]>=options['2pt']['split_val'])
		print('%3.3f percent split'%(data[name][mask].size*100./data[name].size))
	else:
		print('No catalogue split required.')
		mask = np.ones(data.size).astype(bool)

	print('Setting up correlations')
	cat1 = treecorr.Catalog(x=data['x'][mask], y=data['y'][mask])
	
	if splitflag:
		cat2 = treecorr.Catalog(x=data['x'][np.invert(mask)], y=data['y'][np.invert(mask)])

		c1c1 = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=binning)
		c2c2 = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=binning)
		c1c2 = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=binning)
		c2c1 = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=binning)
		suffix = '_splitby%s'%options['2pt']['split']

		rcat11, rcat12, rcat21, rcat22 = randoms(cat1, cat2, period=periods[options['simulation']])

	else:
		rcat11, rcat12, rcat21, rcat22 = randoms(cat1, cat1, period=periods[options['simulation']])
		suffix=''

	if (nlow>=0):
		suffix+='-ndm_part_low%d'%nlow
	if (nhigh>0):
		suffix+='-ndm_part_high%d'%nhigh
	

	print('Computing correlation functions.')

	if splitflag:
		print('11')
		c1c1.process(cat1,cat1)
		c1c1 = finish_nn(c1c1, cat1, cat1, rcat11, rcat12, options, binning)
		print('22')
		c2c2.process(cat2,cat2)
		c2c2 = finish_nn(c2c2, cat2, cat2, rcat21, rcat22, options, binning)	
		print('12')
		c1c2.process(cat1,cat2)
		c1c2 = finish_nn(c1c2, cat1, cat2, rcat11, rcat21, options, binning)
		print('12')
		c2c1.process(cat2,cat1)
		c2c1 = finish_nn(c2c1, cat2, cat1, rcat21, rcat11, options, binning)

		if options['2pt']['errors']:
			dc1c1 = errors.jackknife(data[mask], data[mask], options, nbins=binning)
			dc2c2 = errors.jackknife(data[np.invert(mask)], data[np.invert(mask)], options, nbins=binning)
			dc1c2 = errors.jackknife(data[mask], data[np.invert(mask)], options, nbins=binning)
			dc2c1 = errors.jackknife(data[np.invert(mask)], data[mask], options, nbins=binning)

		else:
			dc1c1 = np.zeros(c1c1.xi.size)
			dc2c2 = np.zeros(c2c2.xi.size)
			dc1c2 = np.zeros(c1c2.xi.size)
			dc2c1 = np.zeros(c2c1.xi.size)

		export_treecorr_output('%s/gg_proj_corr_11%s.txt'%(options['2pt']['savedir'], suffix), c1c1, dc1c1)
		export_treecorr_output('%s/gg_proj_corr_22%s.txt'%(options['2pt']['savedir'], suffix), c2c2, dc2c2)
		export_treecorr_output('%s/gg_proj_corr_12%s.txt'%(options['2pt']['savedir'], suffix), c1c2, dc1c2)
		export_treecorr_output('%s/gg_proj_corr_21%s.txt'%(options['2pt']['savedir'], suffix), c2c1, dc2c1)

	cat0 = treecorr.Catalog(x=data['x'], y=data['y'])
	c0c0 = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=binning)

	rcat01, rcat02, rcat21, rcat22 = randoms(cat0, cat0, period=periods[options['simulation']])

	print('00')
	c0c0.process(cat0,cat0)
	c0c0 = finish_nn(c0c0, cat0, cat0, rcat01, rcat02, options, binning)
	if options['2pt']['errors']:
		dc0c0 = errors.jackknife(data, data, options, nbins=binning)
	else:
		dc0c0 = np.zeros(c0c0.xi.size)

	export_treecorr_output('%s/gg_proj_corr_00%s.txt'%(options['2pt']['savedir'], suffix), c0c0, dc0c0)	

	print('Done')

def randoms(cat1, cat2, period=100):
	# Initialise randoms
	# Fix the random seed so the randoms are always the same for a catalog of given length
	# Should probably migrate this to a config option
	np.random.seed(9000)
	rcat11 = treecorr.Catalog(x=np.random.rand(cat1.ntot)*period, y=np.random.rand(cat1.ntot)*period)
	rcat12 = treecorr.Catalog(x=np.random.rand(cat1.ntot)*period, y=np.random.rand(cat1.ntot)*period)
	rcat21 = treecorr.Catalog(x=np.random.rand(cat2.ntot)*period, y=np.random.rand(cat2.ntot)*period)
	rcat22 = treecorr.Catalog(x=np.random.rand(cat2.ntot)*period, y=np.random.rand(cat2.ntot)*period)

	return rcat11, rcat12, rcat21, rcat22

def randoms_halotools(cat1, period=100):
	# Initialise randoms
	# Fix the random seed so the randoms are always the same for a catalog of given length
	# Should probably migrate this to a config option
	np.random.seed(9000)
	rcat = np.array([np.random.rand(cat1.size)*period, np.random.rand(cat1.size)*period, np.random.rand(cat1.size)*period])
	return rcat.T

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

def export_treecorr_output(filename,corr,errors):
    R = np.exp(corr.logr)

    xi = corr.xi

    out = np.vstack((R,xi,corr.weight, errors))
    print('Saving %s'%filename)
    np.savetxt(filename, out.T)

def compute(options, binning):
	if 'positions' in options['2pt'].keys():
		posfile = options['2pt']['positions']
		data = fi.FITS(posfile)[-1].read()
	else:
		posfile = options['2pt']['shapes']
		data = fi.FITS(posfile)[-1].read()

	print('Position data : %s'%posfile)

	splitflag=options['2pt']['split']
	if 'npart_cut' in options['2pt'].keys():
		nlow = options['2pt']['npart_cut']
		print('Imposing additional cut at npart>%d'%nlow)
		data = data[(data['npart_dm']>nlow)]
	else:
		nlow=-1

	if 'npart_cut_upper' in options['2pt'].keys():
		nhigh = options['2pt']['npart_cut_upper']
		print('Imposing additional cut at npart<%d'%nhigh)
		data = data[(data['npart_dm']<nhigh)]
	else:
		nhigh = -1

	if splitflag is not None:
		name = options['2pt']['split']
		print('Dividing catalogue by %s'%name)
		mask = (data[name]>=options['2pt']['split_val'])
		print('%3.3f percent split'%(data[name][mask].size*100./data[name].size))
	else:
		print('No catalogue split required.')
		mask = np.ones(data.size).astype(bool)

	print('Setting up correlations')
	cat1 = data[mask]
	rbins = np.logspace(np.log10(options['2pt']['rmin']), np.log10(options['2pt']['rmax']), binning+1)

	if splitflag:
		cat2 = data[np.invert(mask)]
		suffix = '_splitby%s'%options['2pt']['split']
	else:
		suffix=''

	if (nlow>=0):
		suffix+='-ndm_part_low%d'%nlow
	if (nhigh>0):
		suffix+='-ndm_part_high%d'%nhigh
	

	print('Computing correlation functions.')

	if splitflag:
		print('11')
		c1c1 = compute_gg(cat1,cat1,options,binning)
		print('22')
		c2c2 = compute_gg(cat2,cat2,options, binning)
		print('12')
		c1c2 = compute_gg(cat1,cat2,options, binning)
		print('12')
		c2c1 = compute_gg(cat2,cat1,options, binning)

		if options['2pt']['errors']:
			dc1c1 = errors.jackknife(data[mask], data[mask], options)
			dc2c2 = errors.jackknife(data[np.invert(mask)], data[np.invert(mask)], options)
			dc1c2 = errors.jackknife(data[mask], data[np.invert(mask)], options)
			dc2c1 = errors.jackknife(data[np.invert(mask)], data[mask], options)

		else:
			dc1c1 = np.zeros(binning)
			dc2c2 = np.zeros(binning)
			dc1c2 = np.zeros(binning)
			dc2c1 = np.zeros(binning)

		export('%s/gg_proj_corr_11%s-ht.txt'%(options['2pt']['savedir'], suffix), rbins, c1c1, dc1c1)
		export('%s/gg_proj_corr_22%s-ht.txt'%(options['2pt']['savedir'], suffix), rbins, c2c2, dc2c2)
		export('%s/gg_proj_corr_12%s-ht.txt'%(options['2pt']['savedir'], suffix), rbins, c1c2, dc1c2)
		export('%s/gg_proj_corr_21%s-ht.txt'%(options['2pt']['savedir'], suffix), rbins, c2c1, dc2c1)

	cat0 = data
	c0c0 = compute_gg(cat0, cat0, options, binning)

	print('00')
	if options['2pt']['errors']:
		dc0c0 = errors.jackknife(data, data, options)
	else:
		dc0c0 = np.zeros(binning)

	export('%s/gg_proj_corr_00%s-ht.txt'%(options['2pt']['savedir'], suffix), rbins, c0c0, dc0c0)
		

	print('Done')


def compute_gg(cat1, cat2, options, nbins):
	pvec1 = np.vstack((cat1['x'], cat1['y'],  cat1['z'])).T
	pvec2 = np.vstack((cat2['x'], cat2['y'],  cat2['z'])).T
	rbins = np.logspace(np.log10(options['2pt']['rmin']), np.log10(options['2pt']['rmax']), nbins+1 )
	rvec2 = randoms_halotools(cat2)

	pi_max = options['2pt']['pi_max'] 

	gg = wp(pvec2, rbins, pi_max, sample2=pvec1, randoms=rvec2, num_threads=1, estimator='Landy-Szalay') 
	#period=periods[options['simulation']]

	return gg


def export(path, edges, result, error):
	x = np.sqrt(edges[:-1]*edges[1:])
	out = np.vstack((x, result, error))
	print('Exporting to %s'%path)
	np.savetxt(path, out.T)