import fitsio as fi
import numpy as np
import glob
import copy
import argparse
import os 
import pylab as plt
plt.switch_backend('agg')

parser = argparse.ArgumentParser()
parser.add_argument('-f','--filename', action='store', type=str)
parser.add_argument('--nofits', action='store_true')
parser.add_argument('--nofig', action='store_true')
args = parser.parse_args()

files = glob.glob(args.filename)
print 'found %d files to process'%len(files)

baryons1=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat-nthreshold5.fits')['baryons'][:]
dm1=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat-nthreshold5.fits')['dm'][:]

select = (dm1['npart']>1000) & (baryons1['npart']>0) & (np.isfinite(baryons1['x']) & np.isfinite(baryons1['y']) & np.isfinite(baryons1['z']))

dat = copy.deepcopy(baryons1[select])
for f in files:
	d=fi.FITS(f)['baryons'].read()
	for comp in ['x', 'y', 'z']: dat[comp][(d[comp]!=baryons1[comp][select])] = d[comp][(d[comp]!=baryons1[comp][select])]
	print f

symmetrised=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat-nthreshold5.fits')['baryons'][:]
symmetrised['x'][select]=dat['x']
symmetrised['y'][select]=dat['y']
symmetrised['z'][select]=dat['z']

import pdb ; pdb.set_trace()
assert (symmetrised['x'][select]!=baryons1['x'][select]).all()


if not args.nofits:
	os.system('rm symmetrised_catalogue.fits')
	outfits = fi.FITS('symmetrised_catalogue.fits', 'rw')
	outfits.write(symmetrised)
	outfits.close()


if args.nofig:
	exit()

from mbii.readsubhalo import *
from mbii.properties import *
root_folder='/physics/yfeng1/mb2'
snapshot='085'

snap = SnapDir(snapshot, root_folder)
h = snap.readsubhalo()

colours=['purple', 'pink', 'plum', 'hotpink', 'forestgreen', 'orange', 'k', 'gray']*100
group_ids = np.array(h['groupid'])[select]
groups = np.unique(group_ids)

nhalo=200
print 'will map out the most massive %d halos'%nhalo

for i, g in enumerate(groups):
	if i>nhalo: break
	plt.plot(symmetrised['x'][select][(group_ids==g)], symmetrised['y'][select][(group_ids==g)], '.', color=colours[i])

plt.xlim(0,100000) ; plt.ylim(0,100000) ; plt.savefig('/home/ssamurof/symmetrised_halos.png')