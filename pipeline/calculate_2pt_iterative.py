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
	print('Using MPI (%d processes)'%size)
else:
	rank = 0
	size = 1


options = yaml.load(open(args.config))
binning = lib.parse_binning(options)

# Fixed values from the snapshot 85 catalogue
# should add approximately equal numbers of extra objects for each step  
Nbins = 10**np.array([3.27783833, 3.45742769, 3.65705585, 3.98793427])
	#300,645,1239,1954,2700,3775])
correlations = options['2pt']['ctypes'].split()
options['2pt']['split'] = None

mode = options['2pt']['mode']

for correlation in correlations:
	print('Processing %s'%correlation )
	if (mode=='errors'):
		exec('from mbii.pipeline.twopoint import calculate_%s_errors as fns'%correlation)
	else:
		exec('from mbii.pipeline.twopoint import calculate_%s as fns'%correlation)
	for nlow in Nbins:
		options['2pt']['npart_cut'] = nlow
		nbins = binning[correlation]
		fns.compute(options, nbins)