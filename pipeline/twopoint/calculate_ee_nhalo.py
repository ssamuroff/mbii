import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
#import mbii.lego_tools as util
from halotools.mock_observables.alignments import ee_3d_one_two_halo_decomp as ee_3d
from mbii.pipeline.twopoint.jackknife import ee_nhalo as errors


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
	c1c1_1h, c1c1_2h = compute_ee(cat1, cat1, options) 

	if options['2pt']['errors']:
		dc1c1_1h, dc1c1_2h = errors.jackknife(data[mask], data[mask], options, rbins=rbins)
	else:
		dc1c1_1h, dc1c1_2h = np.zeros(c1c1.size), np.zeros(c1c1.size)

	if splitflag:
		cat2 = data[np.invert(mask)]
		suffix = '_splitby%s'%options['2pt']['split']
	else:
		suffix=''

	if splitflag:
		print '21'
		c2c1_1h, c2c1_2h = compute_ee(cat2, cat1, options, rbins=rbins)

		print '22'
		c2c2_1h, c2c2_2h = compute_ee(cat2, cat2, options, rbins=rbins)
	
		print '12'
		c1c2_1h, c1c2_2h = compute_ee(cat1, cat2, options, rbins=rbins)
		

		if options['2pt']['errors']:
			dc2c2_1h, dc2c2_2h = errors.jackknife(cat2, cat2, options, rbins=rbins)
			dc1c2_1h, dc1c2_2h = errors.jackknife(cat1, cat2, options, rbins=rbins)
			dc2c1_1h, dc2c1_2h = errors.jackknife(cat2, cat1, options, rbins=rbins)

		else:
			dc1c1_1h, dc1c1_2h = np.zeros(c1c1.size), np.zeros(c1c1.size)
			dc2c2_1h, dc2c2_2h = np.zeros(c2c2.size), np.zeros(c2c2.size)
			dc1c2_1h, dc1c2_2h = np.zeros(c1c2.size), np.zeros(c1c2.size)
			dc2c1_1h, dc2c1_2h = np.zeros(c2c1.size), np.zeros(c2c1.size)

		
		export_array('%s/EE_corr_11%s_1h.txt'%(options['2pt']['savedir'], suffix), rbins, c1c1_1h, dc1c1_1h)
		export_array('%s/EE_corr_22%s_1h.txt'%(options['2pt']['savedir'], suffix), rbins, c2c2_1h, dc2c2_1h)
		export_array('%s/EE_corr_12%s_1h.txt'%(options['2pt']['savedir'], suffix), rbins, c1c2_1h, dc1c2_1h)
		export_array('%s/EE_corr_21%s_1h.txt'%(options['2pt']['savedir'], suffix), rbins, c2c1_1h, dc2c1_1h)

		export_array('%s/EE_corr_11%s_2h.txt'%(options['2pt']['savedir'], suffix), rbins, c1c1_2h, dc1c1_2h)
		export_array('%s/EE_corr_22%s_2h.txt'%(options['2pt']['savedir'], suffix), rbins, c2c2_2h, dc2c2_2h)
		export_array('%s/EE_corr_12%s_2h.txt'%(options['2pt']['savedir'], suffix), rbins, c1c2_2h, dc1c2_2h)
		export_array('%s/EE_corr_21%s_2h.txt'%(options['2pt']['savedir'], suffix), rbins, c2c1_2h, dc2c1_2h)

	print '00'
	cat0 = data
	c0c0_1h, c0c0_2h = compute_ee(cat0, cat0, options, rbins=rbins)
	
	if options['2pt']['errors']:
		dc0c0_1h, dc0c0_2h = errors.jackknife(data, data, options, rbins=rbins)
	else:
		dc0c0_1h, dc0c0_2h = np.zeros(c0c0.xi.size), np.zeros(c0c0.xi.size)
	export_array('%s/EE_corr_00%s_1h.txt'%(options['2pt']['savedir'], suffix), rbins, c0c0_1h, dc0c0_1h)
	export_array('%s/EE_corr_00%s_2h.txt'%(options['2pt']['savedir'], suffix), rbins, c0c0_2h, dc0c0_2h)
		

	print 'Done'

def compute_ee(cat1, cat2, options, period=100., rbins=None):
	avec2 = np.vstack((cat2['a1'], cat2['a2'], cat2['a3'])).T
	avec1 = np.vstack((cat1['a1'], cat1['a2'], cat1['a3'])).T
	pvec1 = np.vstack((cat1['x'], cat1['y'], cat1['z'])).T
	pvec2 = np.vstack((cat2['x'], cat2['y'], cat2['z'])).T
	hids1 = cat1['halo_id']
	hids2 = cat2['halo_id']

	if options['2pt']['binning']=='log' and (rbins is None):
		rbins = np.logspace(np.log10(options['2pt']['rmin']), np.log10(options['2pt']['rmax']), options['2pt']['nbin'])
	elif options['2pt']['binning']=='equal' and (rbins is None):
		rbins = util.equalise_binning(cat1,cat2,options['2pt']['rmin'], options['2pt']['rmax'], options['2pt']['nbin'])

	ee_1h, ee_2h = ee_3d(pvec2, avec2, hids2, pvec1, avec1, hids1, rbins, period=period, num_threads=1) 

	return ee_1h, ee_2h

def export_array(path, edges, result, error):
	x = np.sqrt(edges[:-1]*edges[1:])
	out = np.vstack((x, result, error))
	print 'Exporting to %s'%path
	np.savetxt(path, out.T)

