import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
#import mbii.lego_tools as util
from mbii.pipeline.twopoint.jackknife import gg as errors 

def generate_random_particle_cat(options, npart=2911409):
    print('Generating random subsample of %d dark matter particles'%npart)
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
	particles = generate_random_particle_cat(options)

	splitflag=options['2pt']['split']
	nlow=-1
	nhigh = -1

	print('No catalogue split required.')
	mask = np.ones(particles['x'].size).astype(bool)

	print('Setting up correlations')
	cat1 = treecorr.Catalog(x=particles['x'][mask], y=particles['y'][mask], z=particles['z'][mask])
	c1c1 = treecorr.NNCorrelation(min_sep=options['2pt']['rmin'], max_sep=options['2pt']['rmax'], nbins=options['2pt']['nbin'])

	rcat11, rcat12, rcat21, rcat22 = randoms(cat1, cat1)
	suffix=''

	print('Computing correlation functions.')
	print('11')
	c1c1.process(cat1,cat1)
	c1c1 = finish_nn(c1c1, cat1, cat1, rcat11, rcat12, options)
	if options['2pt']['errors']:
		dc1c1 = errors.jackknife(particles[mask], particles[mask], options)
	else:
		dc1c1 = np.zeros(c1c1.xi.size)

	export_treecorr_output('%s/gg_corr_00%s.txt'%(options['2pt']['savedir'], suffix), c1c1, dc1c1)
		

	print('Done')

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
