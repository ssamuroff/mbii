import numpy as np
import matplotlib.pyplot as plt
from mbii.whizzy_plot import *
plt.switch_backend('agg')
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d
from mbii.mbdb import mbdb
import mbii.lego_tools as util
import mbii
import treecorr
import fitsio as fi
import halotools
import time, argparse

from mbii.readsubhalo import *
from mbii.properties import *


parser = argparse.ArgumentParser()
parser.add_argument('--mpi', action='store_true')
args = parser.parse_args()

if args.mpi:
	import mpi4py.MPI
	rank = mpi4py.MPI.COMM_WORLD.Get_rank()
	size = mpi4py.MPI.COMM_WORLD.Get_size()
	print 'Using MPI (%d processes)'%size
else:
	rank = 0
	size = 1

data = fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat-nthreshold5.fits')['baryons'][:]
dm = fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat-nthreshold5.fits')['dm'][:]

root_folder='/physics/yfeng1/mb2'
snapshot='085'

snap = SnapDir(snapshot, root_folder)
h = snap.readsubhalo()
select = (dm['npart']>1000) & (data['npart']>0) & (np.isfinite(h['pos'].T[0]) & np.isfinite(h['pos'].T[1]) & np.isfinite(h['pos'].T[2]))
flags = mbii.mbdb.identify_centrals(h, mask=select, rank=rank, size=size)
