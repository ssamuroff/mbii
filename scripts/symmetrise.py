import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from mbii.mbdb import mbdb
import mbii
import treecorr
import fitsio as fi
import sys

import mbii.lego_tools as utils
import  mbii.properties as prop
import argparse

root_folder='/physics/yfeng1/mb2'


parser = argparse.ArgumentParser()
parser.add_argument('-c','--catalogue', action='store', type=str)
parser.add_argument('--component', action='store', type=str)
parser.add_argument('--nrot', action='store', type=int)
parser.add_argument('--snapshot', action='store', default='085', type=int)
args = parser.parse_args()

snapshot=args.snapshot
print 'Will process snapshot %s'%snapshot
print 'Matter component:', args.component
print 'Input catalogue:', args.catalogue
print 'Using %d random rotations per halo'%args.nrot

data = fi.FITS(args.catalogue)[args.component][:]

dm = fi.FITS(args.catalogue)['dm'][:]

select = (dm['npart']>1000) & (data['npart']>0) & (np.isfinite(data['x']) & np.isfinite(data['y']) & np.isfinite(data['z']))
utils.symmetrise_catalogue(data, mask=select, nrot=args.nrot)
