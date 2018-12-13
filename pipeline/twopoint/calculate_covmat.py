import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
#import mbii.lego_tools as util
from mbii.pipeline.twopoint.jackknife import covmat as errors 

period={'massiveblackii':100, 'illustris':75}

def compute(options, binning, snapshots, mpi=False):

	if mpi:
		import mpi4py.MPI
		rank = mpi4py.MPI.COMM_WORLD.Get_rank()
		nthread = mpi4py.MPI.COMM_WORLD.Get_size()
	else:
		rank = 0
		nthread = 1

	data = load_catalogues(options, snapshots)

	dc0c0 = errors.jackknife(options['covariance']['ctypes'].split(), data, data, options, nbins=binning, rank=rank, nthread=nthread)
	#np.savetxt('bootstrap-covmat_all.txt',dc0c0)

	print('Done')


def load_catalogues(options, snapshots):

	#print('Shape data : %s'%options['2pt']['shapes'])

	filetemp = options['2pt']['shapes'].split('snapshot')[0]
	data = {}

	for s in snapshots:
		filename = filetemp + 'snapshot%d.fits'%int(s)
		print('Snapshot %d: %s'%(int(s), filename))
		data[int(s)] = fi.FITS(filename)[-1].read()

	return data

def export_array(path, edges, result, error):
	x = np.sqrt(edges[:-1]*edges[1:])
	out = np.vstack((x, result, error))
	print('Exporting to %s'%path)
	np.savetxt(path, out.T)

