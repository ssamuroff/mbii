import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
#import mbii.lego_tools as util
from halotools.mock_observables.alignments import gi_plus_projected
from mbii.pipeline.twopoint.jackknife import giplus_proj as errors 

period={'massiveblackii':100, 'illustris':75}

def compute(options, binning):

	shapefile = options['2pt']['shapes']

	if 'positions' in options['2pt'].keys():
		posfile = options['2pt']['positions']
	else:
		posfile = shapefile

	print('Shape data : %s'%shapefile)
	print('Density data : %s'%posfile)

	shapes = fi.FITS(shapefile)[-1].read()
	positions = fi.FITS(posfile)[-1].read()

	if 'npart_cut' in options['2pt'].keys():
		nlow = options['2pt']['npart_cut']
		print('Imposing additional cut at npart>%d'%nlow)
		shapes = shapes[(shapes['npart_dm']>nlow)]
		if (shapefile==posfile):
			positions = positions[positions['npart_dm']>nlow]
	else:
		nlow = -1

	if 'npart_cut_upper' in options['2pt'].keys():
		nhigh = options['2pt']['npart_cut_upper']
		print('Imposing additional cut at npart<%d'%nhigh)
		shapes = shapes[(shapes['npart_dm']<nhigh)]
		if (shapefile==posfile):
			positions = positions[positions['npart_dm']<nhigh]
	else:
		nhigh = -1

	splitflag=options['2pt']['split']

	if splitflag is not None:
		name = options['2pt']['split']
		print('Dividing catalogue by %s'%name)
		mask = (shapes[name]>=options['2pt']['split_val'])
		print('%3.3f percent split'%(shapes[name][mask].size*100./shapes[name].size))
	else:
		print('No catalogue split required.')
		mask = np.ones(shapes.size).astype(bool)

	rpbins = np.logspace(np.log10(options['2pt']['rmin']), np.log10(options['2pt']['rmax']), binning+1)

	print('Setting up correlations')

	if splitflag:
		suffix = '_splitby%s'%options['2pt']['split']
	else:
		suffix=''

	if (nlow>=0):
		suffix+='-ndm_part_low%d'%nlow
	if (nhigh>0):
		suffix+='-ndm_part_high%d'%nhigh

	# don't need randoms here if we know the period of the box
	print('Computing correlation functions.')

	if splitflag:
		cat1 = shapes[mask]
		cat2 = shapes[np.invert(mask)]

		print('11')
		c1c1 = compute_giplus(cat1, cat1, options, period=period[options['simulation']], nbins=binning)

		print('22')
		c2c2 = compute_giplus(cat2, cat2, options, rpbins=rpbins, period=period[options['simulation']], nbins=binning)
	
		print('12')
		c1c2 = compute_giplus(cat1, cat2, options, rpbins=rpbins, period=period[options['simulation']], nbins=binning)
		print('21')
		c2c1 = compute_giplus(cat2, cat1, options, rpbins=rpbins, period=period[options['simulation']], nbins=binning)

		if options['2pt']['errors']:
			dc1c1 = errors.jackknife(data[mask], data[mask], options, nbins=binning)
			dc2c2 = errors.jackknife(cat2, cat2, options, nbins=binning)
			dc1c2 = errors.jackknife(cat1, cat2, options, nbins=binning)
			dc2c1 = errors.jackknife(cat2, cat1, options, nbins=binning)

		else:
			dc1c1 = np.zeros(c1c1.size)
			dc2c2 = np.zeros(c2c2.size)
			dc1c2 = np.zeros(c1c2.size)
			dc2c1 = np.zeros(c2c1.size)

		export_array('%s/GIplus_proj_corr_11%s.txt'%(options['2pt']['savedir'], suffix), rpbins, c1c1, dc1c1) 
		export_array('%s/GIplus_proj_corr_22%s.txt'%(options['2pt']['savedir'], suffix), rpbins, c2c2, dc2c2)
		export_array('%s/GIplus_proj_corr_12%s.txt'%(options['2pt']['savedir'], suffix), rpbins, c1c2, dc1c2)
		export_array('%s/GIplus_proj_corr_21%s.txt'%(options['2pt']['savedir'], suffix), rpbins, c2c1, dc2c1)

	print('00')
	c0c0 = compute_giplus(positions, shapes, options, rpbins=rpbins, period=period[options['simulation']], nbins=binning)
	
	if options['2pt']['errors']:
		dc0c0 = errors.jackknife(positions, shapes, options, nbins=binning)
	else:
		dc0c0 = np.zeros(binning)
	export_array('%s/GIplus_proj_corr_00%s.txt'%(options['2pt']['savedir'], suffix), rpbins, c0c0, dc0c0)
		

	print('Done')

def compute_giplus(cat1, cat2, options, period=100., rpbins=None, nbins=6):

	aname = 'a%d'
	ename = 'e%d'
	if ('shapes_suffix' in options['2pt'].keys()):
		aname+=options['2pt']['shapes_suffix']
		ename+=options['2pt']['shapes_suffix']

	avec = np.vstack((cat2[aname%1], cat2[aname%2])).T
	pvec1 = np.vstack((cat1['x'], cat1['y'], cat1['z'])).T
	pvec2 = np.vstack((cat2['x'], cat2['y'], cat2['z'])).T
	evec = np.sqrt(cat2[ename%1]*cat2[ename%1] + cat2[ename%2]*cat2[ename%2])

	rpbins = np.logspace(np.log10(options['2pt']['rmin']), np.log10(options['2pt']['rmax']), nbins+1)
	
	pi_max = options['2pt']['pi_max']

	mask2=avec.T[0]!=0.0

	gip = gi_plus_projected(pvec2[mask2], avec[mask2], evec[mask2], pvec1, rpbins, pi_max, period=period, num_threads=1) 

	return gip

def export_array(path, edges, result, error):
	x = np.sqrt(edges[:-1]*edges[1:])
	out = np.vstack((x, result, error))
	print('Exporting to %s'%path)
	np.savetxt(path, out.T)

