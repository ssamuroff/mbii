import argparse
import numpy as np
import yaml
import fitsio as fi
from mbii.pipeline.twopoint.jackknife import differr 


parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--config', '-c', type=str, action='store')
args = parser.parse_args()

options = yaml.load(open(args.config))


# no longer care about making this nice to use
# doing something that works in as little time as possible
fid = fi.FITS('/Users/hattifattener/Documents/ias/mbii/cats/final2/postprocessed/base_subhalo_shapes-v10-ndm1000-nst300-snapshot85.fits')[-1].read()
sym = fi.FITS('/Users/hattifattener/Documents/ias/mbii/cats/final2/postprocessed/base_subhalo_shapes-v10-ndm1000-nst300-snapshot85-symmetrised0.fits')[-1].read()

# Three permutations of the symmetrised catalogue:
# Using an alternative definition for central flagging
p1 = fid
p1_sym = fi.FITS('/Users/hattifattener/Documents/ias/mbii/cats/final2/postprocessed/base_subhalo_shapes-v10-ndm1000-nst300-snapshot85-most_massive_cent_symmetrised0.fits')[-1].read()

# Rotating about the central galaxy, not the halo centroid
p2 = fid
p2_sym = fi.FITS('/Users/hattifattener/Documents/ias/mbii/cats/final2/postprocessed/base_subhalo_shapes-v10-ndm1000-nst300-snapshot85-galaxy_symmetrised0.fits')[-1].read()

# Using dark matter shapes rather than the baryonic component
p3 = fi.FITS('/Users/hattifattener/Documents/ias/mbii/cats/final2/postprocessed/base_subhalo_shapes-v10-ndm1000-nst300-snapshot85.fits')[-1].read()
p3_sym = fi.FITS('/Users/hattifattener/Documents/ias/mbii/cats/final2/postprocessed/base_subhalo_shapes-v10-ndm1000-nst300-snapshot85-symmetrised0.fits')[-1].read()

for name in ['e1','e2','a1','a2','a3']:
	p3[name]=p3['%s_dm'%name]
	p3_sym[name]=p3_sym['%s_dm'%name]

correlations = ['ed', 'ee', 'gi_plus_projected', 'ii_plus_projected']
for correlation in correlations:
	print('Correlation : %s'%correlation)

	print('Processing permutation 1.')
	r, df, vardf = differr.jackknife(correlation, [fid,p1], [fid,p1], [sym,p1_sym], [sym,p1_sym], options)
	differr.export_array('%s_corrvar_altcent.txt'%correlation, r, df, vardf)

	print('Processing permutation 2.')
	r, df, vardf = differr.jackknife(correlation, [fid,p2], [fid,p2], [sym,p2_sym], [sym,p2_sym], options)
	differr.export_array('%s_corrvar_galaxygalaxysym.txt'%correlation, r, df, vardf)

	print('Processing permutation 3.')
	r, df, vardf = differr.jackknife(correlation, [fid,p3], [fid,p3], [sym,p3_sym], [sym,p3_sym], options)
	differr.export_array('%s_corrvar_dmshapes.txt'%correlation, r, df, vardf)

print('Done.')

