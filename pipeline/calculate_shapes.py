import numpy as np
import yaml
import argparse
from mbii.mbdb import mbdb
import mbii.lego_tools as util
import mbii
import fitsio as fi
import time
import sys

from mbii.readsubhalo import *
from mbii.properties import *

import mbii.shapes_lib as lib


parser = argparse.ArgumentParser()
parser.add_argument('--config', action='store', type=str)
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

options = yaml.load(open(args.config))
simulation = options['simulation']
snapshot = options['catalogues']['snapshot']
root_folder = options['root_folder']

print 'Using simulation : %s' %simulation
print 'Using snapshot : %d'%snapshot

# We have a few possible permutations for the subhalo shapes here
# (a) 2D or 3D shapes? It's the same method really, but the former involves cutting out a 2x2 corner of the 3x3 inertia tensor matrix.
# (b) simple or reduced inertia tensor? The difference is a radial weighting term, which should reduce the impact of stuff on the edges of the profile.
# (c) Inertia tensor or spin?

if (options['catalogues']['shapes_method']=='inertia_tensor'):
    lib.compute_inertia_tensors(options, rank, size, reduced=False, snapshot=snapshot, savedir=options['catalogues']['shapes_dir'])
elif (options['catalogues']['shapes_method']=='reduced_inertia_tensor'):
    lib.compute_inertia_tensors(options, rank, size, reduced=True, snapshot=snapshot, savedir=options['catalogues']['shapes_dir'])
elif (options['catalogues']['shapes_method']=='spin'):
    lib.compute_spin(options, rank, size, component='baryons', snapshot=snapshot, savedir=options['catalogues']['shapes_dir'])
