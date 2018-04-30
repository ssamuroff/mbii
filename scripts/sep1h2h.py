import numpy as np
import fitsio as fi

baryons1=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat-nthreshold5.fits')['baryons'][:]
baryons2=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat_reduced.fits')['baryons'][:]
dm1=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat-nthreshold5.fits')['dm'][:]
dm2=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat_reduced.fits')['dm'][:]
col=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_parent_halo_col.fits')[-1].read()

select = (dm1['npart']>1000) & (baryons1['npart']>0) & (np.isfinite(baryons1['x']) & np.isfinite(baryons1['y']) & np.isfinite(baryons1['z']))

from mbii.readsubhalo import *
from mbii.properties import *


root_folder='/physics/yfeng1/mb2'
snapshot='085'

snap = SnapDir(snapshot, root_folder)
h = snap.readsubhalo()

pids=np.unique(h['groupid'])
ids=np.array(h['groupid'])
subset = np.random.choice(pids, size=pids.size/2, replace=False)
N=len(baryons1)/2

select = np.ones(ids.size).astype(bool)
val=False
for p in pids:
    print p 
    select[ids==p] = val
    val = np.invert(val)

outfits = fi.FITS('/home/ssamurof/subsamp_1h2h.fits', 'rw')
outfits.write(select.astype(int))
outfits.close()
