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
import time
import sys

from mbii.readsubhalo import *
from mbii.properties import *
import  mbii.properties as prop
import  mbii.properties_projected as prop2d
import argparse

root_folder='/physics/yfeng1/mb2'
snapshot='085'

parser = argparse.ArgumentParser()
parser.add_argument('-t','--type', action='store', type=str)
#parser.add_argument('--masswts', action='store_true')

args = parser.parse_args()

#import pdb ; pdb.set_trace()

snap = SnapDir(snapshot, root_folder)
if args.type=='inertia_tensor':
    prop.compute_inertia_tensors(snap, reduced=False)
elif args.type=='inertia_tensor_projected':
    prop2d.compute_inertia_tensors_projected(snap, reduced=False)
elif args.type=='reduced_inertia_tensor_projected':
    prop2d.compute_inertia_tensors_projected(snap, reduced=True)
elif args.type=='reduced_inertia_tensor':
    prop.compute_inertia_tensors(snap, reduced=True)
elif args.type=='spin':
    prop.compute_spin(snap, component='dm', nsubhalo=20000)
