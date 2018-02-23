import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
import mbii.lego_tools as util
from mbii.pipeline.twopoint.jackknife import ed as errors 

def compute(options):
	print 'Position data : %s'%options['output']
	data = fi.FITS(options['2pt']['shapes'])[-1].read()

	splitflag=options['2pt']['split']

	if splitflag is not None:
		name = options['2pt']['split']
		print 'Dividing catalogue by %s'%name
		mask = (data[name]>=options['2pt']['split_val'])
		print '%3.3f percent split'%(data[name][mask].size*100./data[name].size)
	else:
		print 'No catalogue split required.'
		mask = np.ones(data.size).astype(bool)

	print 'Setting up correlations'
	cat1 = treecorr.Catalog(x=data['x'][mask], y=data['y'][mask], z=data['z'][mask])
	c1c1 = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])

	if splitflag:
		cat2 = treecorr.Catalog(x=data['x'][np.invert(mask)], y=data['y'][np.invert(mask)], z=data['z'][mask])
		c2c2 = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])
		c1c2 = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])
		c2c1 = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])
		suffix = '_splitby%s'%options['2pt']['split']
	else:
		suffix=''

	rcat11, rcat12, rcat21, rcat22 = randoms(cat1, cat2)

	print 'Computing correlation functions.'
	print '11'
	c1c1.process(cat1,cat1)
	c1c1 = finish_nn(c1c1, cat1, cat1, rcat11, rcat12, options)
	if options['2pt']['errors']:
		dc1c1 = errors.jackknife(data[mask], data[mask], options)
	else:
		dc1c1 = np.zeros(c1c1.xi.size)

	util.export_treecorr_output('%s/gg_corr_11%s.txt'%(options['2pt']['savedir'], suffix), c1c1, dc1c1)

	if splitflag:
		print '22'
		c2c2.process(cat2,cat2)
		c2c2 = finish_nn(c2c2, cat2, cat2, rcat21, rcat22, options)	
		print '12'
		c1c2.process(cat1,cat2)
		c1c2 = finish_nn(c1c2, cat1, cat2, rcat11, rcat21, options)
		print '12'
		c2c1.process(cat2,cat1)
		c2c1 = finish_nn(c2c1, cat2, cat1, rcat21, rcat11, options)

		if options['2pt']['errors']:
			dc2c2 = errors.jackknife(data[np.invert(mask)], data[np.invert(mask)], options)
			dc1c2 = errors.jackknife(data[mask], data[np.invert(mask)], options)
			dc2c1 = errors.jackknife(data[np.invert(mask)], data[mask], options)

		else:
			dc1c1 = np.zeros(c1c1.xi.size)
			dc2c2 = np.zeros(c2c2.xi.size)
			dc1c2 = np.zeros(c1c2.xi.size)
			dc2c1 = np.zeros(c2c1.xi.size)

		util.export_treecorr_output('%s/gg_corr_22%s.txt'%(options['2pt']['savedir'], suffix), c2c2, dc2c2)
		util.export_treecorr_output('%s/gg_corr_12%s.txt'%(options['2pt']['savedir'], suffix), c1c2, dc1c2)
		util.export_treecorr_output('%s/gg_corr_21%s.txt'%(options['2pt']['savedir'], suffix), c2c1, dc2c1)

	cat0 = treecorr.Catalog(x=data['x'], y=data['y'], z=data['z'])
	c0c0 = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])

	rcat01, rcat02, rcat21, rcat22 = randoms(cat0, cat0)

	print '00'
	c0c0.process(cat0,cat0)
	c0c0 = finish_nn(c0c0, cat0, cat0, rcat01, rcat02, options)
	if options['2pt']['errors']:
		dc0c0 = errors.jackknife(data, data, options)
	else:
		dc0c0 = np.zeros(c0c0.xi.size)

	util.export_treecorr_output('%s/gg_corr_00%s.txt'%(options['2pt']['savedir'], suffix), c0c0, dc0c0)
		

	print 'Done'

def randoms(cat1, cat2):
	# Initialise randoms
	# Fix the random seed so the randoms are always the same for a catalog of given length
	# Should probably migrate this to a config option
	np.random.seed(9000)
	rcat11 = treecorr.Catalog(x=np.random.rand(cat1.ntot)*100, y=np.random.rand(cat1.ntot)*100, z=np.random.rand(cat1.ntot)*100)
	rcat12 = treecorr.Catalog(x=np.random.rand(cat1.ntot)*100, y=np.random.rand(cat1.ntot)*100, z=np.random.rand(cat1.ntot)*100)
	rcat21 = treecorr.Catalog(x=np.random.rand(cat2.ntot)*100, y=np.random.rand(cat2.ntot)*100, z=np.random.rand(cat2.ntot)*100)
	rcat22 = treecorr.Catalog(x=np.random.rand(cat2.ntot)*100, y=np.random.rand(cat2.ntot)*100, z=np.random.rand(cat2.ntot)*100)

	return rcat11, rcat12, rcat21, rcat22

def finish_nn(corr, cat1, cat2, rcat1, rcat2, options):

	if cat2 is None:
		cat2 = cat1

	rr = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])
	rn = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])
	nr = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])

	# Do the pair counting
	print 'Processing randoms',
	print 'RR',
	rr.process(rcat1,rcat2)
	print 'DR',
	nr.process(cat1,rcat2)
	print 'RD',
	rn.process(rcat1,cat2)

	# Finish off
	xi, var = corr.calculateXi(rr, dr=nr, rd=rn)
	setattr(corr, 'xi', xi)

	return corr


