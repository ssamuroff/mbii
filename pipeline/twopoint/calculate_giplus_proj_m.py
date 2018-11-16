import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
import mbii.shapes_lib as slib
#import mbii.lego_tools as util
from halotools.mock_observables.alignments import gi_plus_projected
from mbii.pipeline.twopoint.jackknife import giplus_proj_m as errors 

def generate_random_particle_cat(options, npart=2911409):
    print 'Generating random subsample of %d dark matter particles'%npart
    root_folder = options['root_folder']
    snapshot = '0%d'%options['catalogues']['snapshot']
    h, x, xb = slib.read_subhalo_data(options['simulation'], int(snapshot), root_folder)
    x0 = []
    y0 = []
    z0 = []
    for i in xrange(npart):
        ihalo = np.random.randint(len(x))
        p = x[ihalo]
        isub = np.random.randint(len(p))
        x0.append(p[isub][0]/1e3)
        y0.append(p[isub][1]/1e3)
        z0.append(p[isub][2]/1e3)
    pvec = np.array(npart, dtype=[('x', float), ('y', float), ('z', float)])
    pvec['x'] = x0
    pvec['y'] = y0
    pvec['z'] = z0
    return pvec

def compute(options):

	binning = options['2pt']['binning']

	data = fi.FITS(options['2pt']['shapes'])[-1].read()

	particles = generate_random_particle_cat(options)
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
		rpbins = np.logspace(np.log10(options['2pt']['rpmin']), np.log10(options['2pt']['rpmax']), options['2pt']['nrpbin'])
	elif binning=='equal':
		rpbins = util.equalise_binning(data[mask], data[mask], options['2pt']['rpmin'], options['2pt']['rpmax'], options['2pt']['nbin'])

	print 'Setting up correlations'

	# don't need randoms here if we know the period of the box
	print 'Computing correlation functions.'
	print '11'
	cat1 = data[mask]
	c1c1 = compute_giplus(particles, cat1, options) 

	if options['2pt']['errors']:
		dc1c1 = errors.jackknife(particles, data[mask], options, rpbins=rpbins)
	else:
		dc1c1 = np.zeros(c1c1.size)

	if splitflag:
		cat2 = data[np.invert(mask)]
		suffix = '_splitby%s'%options['2pt']['split']
	else:
		suffix=''

	if splitflag:
		print '22'
		c2c2 = compute_giplus(particles, cat2, options, rpbins=rpbins)
	
		print '12'
		c1c2 = compute_giplus(particles, cat2, options, rpbins=rpbins)
		print '21'
		c2c1 = compute_giplus(particles, cat1, options, rpbins=rpbins)

		if options['2pt']['errors']:
			dc2c2 = errors.jackknife(particles, cat2, options, rpbins=rpbins)
			dc1c2 = errors.jackknife(particles, cat2, options, rpbins=rpbins)
			dc2c1 = errors.jackknife(particles, cat1, options, rpbins=rpbins)

		else:
			dc1c1 = np.zeros(c1c1.size)
			dc2c2 = np.zeros(c2c2.size)
			dc1c2 = np.zeros(c1c2.size)
			dc2c1 = np.zeros(c2c1.size)

		
		export_array('%s/GIplus_proj_m_corr_11%s.txt'%(options['2pt']['savedir'], suffix), rpbins, c1c1, dc1c1)
		export_array('%s/GIplus_proj_m_corr_22%s.txt'%(options['2pt']['savedir'], suffix), rpbins, c2c2, dc2c2)
		export_array('%s/GIplus_proj_m_corr_12%s.txt'%(options['2pt']['savedir'], suffix), rpbins, c1c2, dc1c2)
		export_array('%s/GIplus_proj_m_corr_21%s.txt'%(options['2pt']['savedir'], suffix), rpbins, c2c1, dc2c1)

	print '00'
	cat0 = data
	c0c0 = compute_giplus(particles, cat0, options, rpbins=rpbins)
	
	if options['2pt']['errors']:
		dc0c0 = errors.jackknife(particles, data, options, rpbins=rpbins)
	else:
		dc0c0 = np.zeros(c0c0.xi.size)
	export_array('%s/GIplus_proj_m_corr_00%s.txt'%(options['2pt']['savedir'], suffix), rpbins, c0c0, dc0c0)
		

	print 'Done'

def compute_giplus(cat1, cat2, options, period=100., rpbins=None):
	avec = np.vstack((cat2['a1_dm'], cat2['a2_dm'])).T
	pvec1 = np.vstack((cat1['x'], cat1['y'], cat1['z'])).T
	pvec2 = np.vstack((cat2['x'], cat2['y'], cat2['z'])).T
	evec = np.sqrt(cat2['e1_dm']*cat2['e1_dm'] + cat2['e2_dm']*cat2['e2_dm'])

	if options['2pt']['binning']=='log' and (rpbins is None):
		rpbins = np.logspace(np.log10(options['2pt']['rpmin']), np.log10(options['2pt']['rpmax']), options['2pt']['nrpbin'])
	elif options['2pt']['binning']=='equal' and (rpbins is None):
		rpbins = util.equalise_binning(cat1,cat2,options['2pt']['rpmin'], options['2pt']['rpmax'], options['2pt']['nbin'])
	pi_max = options['2pt']['pi_max']

	gip = gi_plus_projected(pvec2, avec, evec, pvec1, rpbins, pi_max, period=period, num_threads=1) 

	return gip

def export_array(path, edges, result, error):
	x = np.sqrt(edges[:-1]*edges[1:])
	out = np.vstack((x, result, error))
	print 'Exporting to %s'%path
	np.savetxt(path, out.T)

