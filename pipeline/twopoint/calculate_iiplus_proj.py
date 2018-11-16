import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
#import mbii.lego_tools as util
from halotools.mock_observables.alignments import ii_plus_projected
from mbii.pipeline.twopoint.jackknife import iiplus_proj as errors 

period={'massiveblackii':100, 'illustris':75}

def compute(options, binning):
	print('Shape data : %s'%options['2pt']['shapes'])

	data = fi.FITS(options['2pt']['shapes'])[-1].read()
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

	splitflag=options['2pt']['split']

	if splitflag is not None:
		name = options['2pt']['split']
		print('Dividing catalogue by %s'%name)
		mask = (data[name]>=options['2pt']['split_val'])
		print('%3.3f percent split'%(data[name][mask].size*100./data[name].size))
	else:
		print('No catalogue split required.')
		mask = np.ones(data.size).astype(bool)

	if options['2pt']['binning']=='log':
		rpbins = np.logspace(np.log10(options['2pt']['rmin']), np.log10(options['2pt']['rmax']), binning+1)
	elif options['2pt']['binning']=='equal':
		rpbins = util.equalise_binning(data[mask], data[mask], options['2pt']['rmin'], options['2pt']['rmax'], binning+1)

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
		cat1 = data[mask]
		cat2 = data[np.invert(mask)]

		print('11')
		c1c1 = compute_iiplus(cat1, cat1, options, period=period[options['simulation']], nbins=binning)

		print('22')
		c2c2 = compute_iiplus(cat2, cat2, options, rpbins=rpbins, period=period[options['simulation']], nbins=binning)
	
		print('12')
		c1c2 = compute_iiplus(cat1, cat2, options, rpbins=rpbins, period=period[options['simulation']], nbins=binning)
		print('21')
		c2c1 = compute_iiplus(cat2, cat1, options, rpbins=rpbins, period=period[options['simulation']], nbins=binning)

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

		export_array('%s/IIplus_proj_corr_11%s.txt'%(options['2pt']['savedir'], suffix), rpbins, c1c1, dc1c1)
		export_array('%s/IIplus_proj_corr_22%s.txt'%(options['2pt']['savedir'], suffix), rpbins, c2c2, dc2c2)
		export_array('%s/IIplus_proj_corr_12%s.txt'%(options['2pt']['savedir'], suffix), rpbins, c1c2, dc1c2)
		export_array('%s/IIplus_proj_corr_21%s.txt'%(options['2pt']['savedir'], suffix), rpbins, c2c1, dc2c1)

	print('00')
	cat0 = data
	c0c0 = compute_iiplus(cat0, cat0, options, rpbins=rpbins, period=period[options['simulation']], nbins=binning)
	
	if options['2pt']['errors']:
		dc0c0 = errors.jackknife(data, data, options, nbins=binning)
	else:
		dc0c0 = np.zeros(c0c0.xi.size)
	export_array('%s/IIplus_proj_corr_00%s.txt'%(options['2pt']['savedir'], suffix), rpbins, c0c0, dc0c0)	

	print('Done')

def compute_iiplus(cat1, cat2, options, period=100., rpbins=None, nbins=6):

	aname = 'a%d'
	ename = 'e%d'
	if ('shapes_suffix' in options['2pt'].keys()):
		aname+=options['2pt']['shapes_suffix']
		ename+=options['2pt']['shapes_suffix']

	avec1 = np.vstack((cat1[aname%1], cat1[aname%2])).T
	avec2 = np.vstack((cat2[aname%1], cat2[aname%2])).T
	pvec1 = np.vstack((cat1['x'], cat1['y'], cat1['z'])).T
	pvec2 = np.vstack((cat2['x'], cat2['y'], cat2['z'])).T
	evec1 = np.sqrt(cat1[ename%1]*cat1[ename%1] + cat1[ename%2]*cat1[ename%2])
	evec2 = np.sqrt(cat2[ename%1]*cat2[ename%1] + cat2[ename%2]*cat2[ename%2])

	if options['2pt']['binning']=='log' and (rpbins is None):
		rpbins = np.logspace(np.log10(options['2pt']['rmin']), np.log10(options['2pt']['rmax']), nbins+1)
	elif options['2pt']['binning']=='equal' and (rpbins is None):
		rpbins = util.equalise_binning(cat1,cat2,options['2pt']['rmin'], options['2pt']['rmax'], nbins+1)
	pi_max = options['2pt']['pi_max']

	mask1=avec1.T[0]!=0.0
	mask2=avec2.T[0]!=0.0

	iip = ii_plus_projected(pvec2[mask2], avec2[mask2], evec2[mask2], pvec1[mask1], avec1[mask1], evec1[mask1], rpbins, pi_max, period=period, num_threads=1) 

	return iip

def export_array(path, edges, result, error):
	x = np.sqrt(edges[:-1]*edges[1:])
	out = np.vstack((x, result, error))
	print('Exporting to %s'%path)
	np.savetxt(path, out.T)

