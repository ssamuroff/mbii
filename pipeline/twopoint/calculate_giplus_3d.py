import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
#import mbii.lego_tools as util
from halotools.mock_observables.alignments import gi_plus_3d 
from mbii.pipeline.twopoint.jackknife import giplus_3d as errors 

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

	# don't need randoms here if we know the period of the box
	print 'Computing correlation functions.'
	print '11'
	cat1 = data[mask]
	c1c1 = compute_giplus(cat1, cat1, options) 

	if options['2pt']['errors']:
		dc1c1 = errors.jackknife(data[mask], data[mask], options)
	else:
		dc1c1 = np.zeros(c1c1.size)

	if splitflag:
		cat2 = data[np.invert(mask)]
		suffix = '_splitby%s'%options['2pt']['split']
	else:
		suffix=''

	if splitflag:
		print '22'
		c2c2 = compute_giplus(cat2, cat2, options)
	
		print '12'
		c1c2 = compute_giplus(cat1, cat2, options)
		print '21'
		c2c1 = compute_giplus(cat2, cat1, options)

		if options['2pt']['errors']:
			dc2c2 = errors.jackknife(cat2, cat2, options)
			dc1c2 = errors.jackknife(cat1, cat2, options)
			dc2c1 = errors.jackknife(cat2, cat1, options)

		else:
			dc1c1 = np.zeros(c1c1.size)
			dc2c2 = np.zeros(c2c2.size)
			dc1c2 = np.zeros(c1c2.size)
			dc2c1 = np.zeros(c2c1.size)

		rbins = np.logspace(np.log10(options['2pt']['rmin']), np.log10(options['2pt']['rmax']), options['2pt']['nbin'] )
		export_array('%s/GIplus_3d_corr_11%s.txt'%(options['2pt']['savedir'], suffix), rbins, c1c1, dc1c1)
		export_array('%s/GIplus_3d_corr_22%s.txt'%(options['2pt']['savedir'], suffix), rbins, c2c2, dc2c2)
		export_array('%s/GIplus_3d_corr_12%s.txt'%(options['2pt']['savedir'], suffix), rbins, c1c2, dc1c2)
		export_array('%s/GIplus_3d_corr_21%s.txt'%(options['2pt']['savedir'], suffix), rbins, c2c1, dc2c1)

	print '00'
	cat0 = data
	c0c0 = compute_giplus(cat0, cat0, options)
	
	if options['2pt']['errors']:
		dc0c0 = errors.jackknife(data, data, options)
	else:
		dc0c0 = np.zeros(c0c0.xi.size)
	export_array('%s/GIplus_3d_corr_00%s.txt'%(options['2pt']['savedir'], suffix), rbins, c0c0, dc0c0)
		

	print 'Done'

def compute_giplus(cat1, cat2, options):
	avec = np.vstack((cat2['a1'], cat2['a2'], cat2['a3'])).T
	pvec1 = np.vstack((cat1['x'], cat1['y'], cat1['z'])).T
	pvec2 = np.vstack((cat2['x'], cat2['y'], cat2['z'])).T
	evec = np.sqrt(cat2['e1']*cat2['e1'] + cat2['e2']*cat2['e2'])

	rbins = np.logspace(np.log10(options['2pt']['rmin']), np.log10(options['2pt']['rmax']), options['2pt']['nbin'] )

	gip = gi_plus_3d(pvec2, avec, evec, pvec1, rbins, period=100, num_threads=1) 

	return gip


def export_array(path, edges, result, error):
	x = np.sqrt(edges[:-1]*edges[1:])
	out = np.vstack((x, result, error))
	print 'Exporting to %s'%path
	np.savetxt(path, out.T)

