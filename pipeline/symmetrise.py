import numpy as np
import treecorr
import fitsio as fi
import sys, yaml
import mbii.lego_tools as utils
import argparse


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

data = fi.FITS(options['symmetrisation']['catalogue'])[-1].read()

utils.symmetrise_catalogue3(data, seed=options['random_seed'], filename=options['symmetrisation']['catalogue'], savedir=options['symmetrisation']['output'], mask=None, rank=rank, size=size)
