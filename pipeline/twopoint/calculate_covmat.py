import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
#import mbii.lego_tools as util
from mbii.pipeline.twopoint.jackknife import covmat as errors 

period={'massiveblackii':100, 'illustris':75}

def compute(options, binning):
	print('Shape data : %s'%options['2pt']['shapes'])

	data = fi.FITS(options['2pt']['shapes'])[-1].read()

	print('00')
	dc0c0 = errors.jackknife(options['2pt']['ctypes'].split(), data, data, options, nbins=binning)
	np.savetxt('covmat_all.txt',dc0c0)

	print('Done')



def export_array(path, edges, result, error):
	x = np.sqrt(edges[:-1]*edges[1:])
	out = np.vstack((x, result, error))
	print('Exporting to %s'%path)
	np.savetxt(path, out.T)

