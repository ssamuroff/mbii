import numpy as np
import yaml
import fitsio as fi
#import mbii.lego_tools as utils
import mbii.symmetrise_lib as lib
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

nbins = np.array([    1004,     4820,    10919, 46293924])
correlations = options['2pt']['ctypes'].split()
options['2pt']['split'] = None

for correlation in correlations:
	print 'Processing %s'%correlation 
	exec('from mbii.pipeline.twopoint import calculate_%s as fns'%correlation)
	for nlow,nhigh in zip(nbins[:-1], nbins[1:]):
		options['2pt']['npart_cut'] = nlow
		options['2pt']['npart_cut_upper'] = nhigh
		fns.compute(options)