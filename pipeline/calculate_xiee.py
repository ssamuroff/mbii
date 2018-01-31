import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
import mbii.lego_tools as util

def compute(options):
	print 'Shape data : %s'%options['output']
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
	cat1 = treecorr.Catalog(x=data['x'][mask], y=data['y'][mask], z=data['z'][mask], a=data['a1'][mask], b=data['a2'][mask], c=data['a3'][mask])
	c1c1 = treecorr.GGCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])

	if splitflag:
		cat2 = treecorr.Catalog(x=data['x'][np.invert(mask)], y=data['y'][np.invert(mask)], z=data['z'][mask], a=data['a1'][np.invert(mask)], b=data['a2'][np.invert(mask)], c=data['a3'][np.invert(mask)])
		c2c2 = treecorr.GGCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])
		c1c2 = treecorr.GGCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])
		suffix = '_splitby%s'%options['2pt']['split']
	else:
		suffix=''

	print 'Computing correlation functions.'
	c1c1.process(cat1,cat1)
	util.export_treecorr_output('%s/xigg_corr_11%s.txt'%(options['2pt']['savedir'], suffix), c1c1)

	if split:
		c2c2.process(cat2,cat2)
		util.export_treecorr_output('%s/xigg_corr_22%s.txt'%(options['2pt']['savedir'], suffix), c2c2)
		c1c2.process(cat1,cat2)
		util.export_treecorr_output('%s/xigg_corr_12%s.txt'%(options['2pt']['savedir'], suffix), c1c2)

	print 'Done'

