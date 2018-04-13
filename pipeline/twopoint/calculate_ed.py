import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
#import mbii.lego_tools as util
from halotools.mock_observables.alignments import ed_3d
from mbii.pipeline.twopoint.jackknife import ed as errors


def compute(options):
	print 'Shape data : %s'%options['2pt']['shapes']

	binning = options['2pt']['binning']

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

	if binning=='log':
		rbins = np.logspace(np.log10(options['2pt']['rmin']), np.log10(options['2pt']['rmax']), options['2pt']['nbin'])
	elif binning=='equal':
		rbins = util.equalise_binning(data[mask], data[mask], options['2pt']['rmin'], options['2pt']['rmax'], options['2pt']['nbin'])

	print 'Setting up correlations'

	# don't need randoms here if we know the period of the box
	print 'Computing correlation functions.'
	print '11'
	cat1 = data[mask]
	c1c1 = compute_ed(cat1, cat1, options) 

	if options['2pt']['errors']:
		dc1c1 = errors.jackknife(data[mask], data[mask], options, rbins=rbins)
	else:
		dc1c1 = np.zeros(c1c1.size)

	if splitflag:
		cat2 = data[np.invert(mask)]
		suffix = '_splitby%s'%options['2pt']['split']
	else:
		suffix=''

	if splitflag:
		print '22'
		c2c2 = compute_ed(cat2, cat2, options, rbins=rbins)
	
		print '12'
		c1c2 = compute_ed(cat1, cat2, options, rbins=rbins)
		print '21'
		c2c1 = compute_ed(cat2, cat1, options, rbins=rbins)

		if options['2pt']['errors']:
			dc2c2 = errors.jackknife(cat2, cat2, options, rbins=rbins)
			dc1c2 = errors.jackknife(cat1, cat2, options, rbins=rbins)
			dc2c1 = errors.jackknife(cat2, cat1, options, rbins=rbins)

		else:
			dc1c1 = np.zeros(c1c1.size)
			dc2c2 = np.zeros(c2c2.size)
			dc1c2 = np.zeros(c1c2.size)
			dc2c1 = np.zeros(c2c1.size)

		
		export_array('%s/ED_corr_11%s.txt'%(options['2pt']['savedir'], suffix), rbins, c1c1, dc1c1)
		export_array('%s/ED_corr_22%s.txt'%(options['2pt']['savedir'], suffix), rbins, c2c2, dc2c2)
		export_array('%s/ED_corr_12%s.txt'%(options['2pt']['savedir'], suffix), rbins, c1c2, dc1c2)
		export_array('%s/ED_corr_21%s.txt'%(options['2pt']['savedir'], suffix), rbins, c2c1, dc2c1)

	print '00'
	cat0 = data
	c0c0 = compute_ed(cat0, cat0, options, rbins=rbins)
	
	if options['2pt']['errors']:
		dc0c0 = errors.jackknife(data, data, options, rbins=rbins)
	else:
		dc0c0 = np.zeros(c0c0.xi.size)
	export_array('%s/ED_corr_00%s.txt'%(options['2pt']['savedir'], suffix), rbins, c0c0, dc0c0)
		

	print 'Done'

def compute_ed(cat1, cat2, options, period=100., rbins=None):
	avec2 = np.vstack((cat2['a1'], cat2['a2'], cat2['a3'])).T
	avec1 = np.vstack((cat1['a1'], cat1['a2'], cat1['a3'])).T
	pvec1 = np.vstack((cat1['x'], cat1['y'], cat1['z'])).T
	pvec2 = np.vstack((cat2['x'], cat2['y'], cat2['z'])).T

	if options['2pt']['binning']=='log' and (rbins is None):
		rbins = np.logspace(np.log10(options['2pt']['rmin']), np.log10(options['2pt']['rmax']), options['2pt']['nbin'])
	elif options['2pt']['binning']=='equal' and (rbins is None):
		rbins = util.equalise_binning(cat1,cat2,options['2pt']['rmin'], options['2pt']['rmax'], options['2pt']['nbin'])

	ed = ed_3d(pvec2, avec2, pvec1, rbins, period=period, num_threads=1) 

	import pdb ; pdb.set_trace()

	return ed

def export_array(path, edges, result, error):
	x = np.sqrt(edges[:-1]*edges[1:])
	out = np.vstack((x, result, error))
	print 'Exporting to %s'%path
	np.savetxt(path, out.T)

def compute_treecorr(options):
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
	c1c1 = treecorr.NVCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])

	# If required then do jackknife errorbars
	if options['2pt']['errors']:
		dc1c1 = errors.jackknife(data[mask], data[mask], options)
	else:
		dc1c1 = np.zeros(c1c1.xi.size)

	if splitflag:
		cat2 = treecorr.Catalog(x=data['x'][np.invert(mask)], y=data['y'][np.invert(mask)], z=data['z'][mask], a=data['a1'][np.invert(mask)], b=data['a2'][np.invert(mask)], c=data['a3'][np.invert(mask)])
		c2c2 = treecorr.NVCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])
		c1c2 = treecorr.NVCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])
		c2c1 = treecorr.NVCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])
		suffix = '_splitby%s'%options['2pt']['split']
	else:
		suffix=''

	print 'Computing correlation functions.'
	print '11'
	c1c1.process(cat1,cat1)
	export_treecorr_output('%s/ED_corr_11%s.txt'%(options['2pt']['savedir'], suffix), c1c1, dc1c1)

	if splitflag:
		print '22'
		c2c2.process(cat2,cat2)
		if options['2pt']['errors']:
			dc2c2 = errors.jackknife(data[np.invert(mask)], data[np.invert(mask)], options)
			dc1c2 = errors.jackknife(data[mask], data[np.invert(mask)], options)
			dc2c1 = errors.jackknife(data[np.invert(mask)], data[mask], options)

		else:
			dc1c1 = np.zeros(c1c1.xi.size)
			dc2c2 = np.zeros(c2c2.xi.size)
			dc1c2 = np.zeros(c1c2.xi.size)
			dc2c1 = np.zeros(c2c1.xi.size)

		print '12'
		c1c2.process(cat1,cat2)
		print '12'
		c2c1.process(cat2,cat1)

		export_treecorr_output('%s/ED_corr_22%s.txt'%(options['2pt']['savedir'], suffix), c2c2, dc2c2)
		export_treecorr_output('%s/ED_corr_12%s.txt'%(options['2pt']['savedir'], suffix), c1c2, dc1c2)
		export_treecorr_output('%s/ED_corr_21%s.txt'%(options['2pt']['savedir'], suffix), c2c1, dc2c1)

	cat0 = treecorr.Catalog(x=data['x'], y=data['y'], z=data['z'], a=data['a1'], b=data['a2'], c=data['a3'])
	c0c0 = treecorr.NVCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])
	c0c0.process(cat0,cat0)

	if options['2pt']['errors']:
		dc0c0 = errors.jackknife(data, data, options)
	else:
		dc0c0 = np.zeros(c0c0.xi.size)
	export_treecorr_output('%s/ED_corr_00%s.txt'%(options['2pt']['savedir'], suffix), c0c0, dc0c0)

		

	print 'Done'

def export_treecorr_output(filename,corr,errors):
    R = np.exp(corr.logr)

    xi = corr.xi

    out = np.vstack((R,xi,corr.weight, errors))
    print 'Saving %s'%filename
    np.savetxt(filename, out.T)

