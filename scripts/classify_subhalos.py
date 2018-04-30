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

from mbii.readsubhalo import *
from mbii.properties import *


root_folder='/physics/yfeng1/mb2'
snapshot='085'

snap = SnapDir(snapshot, root_folder)
h = snap.readsubhalo()
mbii.mbdb.classify_subhalos(h)
