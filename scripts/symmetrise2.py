import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
#from mbii.mbdb import mbdb
#import mbii
import treecorr
import fitsio as fi
import sys

import mbii.lego_tools as utils
#import  mbii.properties as prop
import argparse

root_folder='/physics/yfeng1/mb2'


parser = argparse.ArgumentParser()
parser.add_argument('-c','--catalogue', action='store', type=str)
parser.add_argument('--component', action='store', type=str)
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

print 'Matter component:', args.component
print 'Input catalogue:', args.catalogue

data = fi.FITS(args.catalogue)[args.component][:]
dm = fi.FITS(args.catalogue)['dm'][:]

#select = np.ones(len(data)).astype(bool)
select = (dm['npart']>1000) & (data['npart']>300) & (np.isfinite(data['x']) & np.isfinite(data['y']) & np.isfinite(data['z'])) & (data['x']<100000) & (data['y']<100000) & (data['z']<100000) & (data['x']>0) & (data['y']>0) & (data['z']>0)
utils.symmetrise_catalogue2(data[select], mask=select, rank=rank, size=size)
