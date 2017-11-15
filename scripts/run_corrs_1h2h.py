import sys
sys.path.append('/home/ssamurof/.local/lib/python2.7/site-packages/')
sys.path.append('/home/ssamurof/.local/lib/python2.6/site-packages/')
import matplotlib.pyplot as plt
from mbii.whizzy_plot import *
plt.switch_backend('agg')
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d
from mbii.mbdb import mbdb
import mbii
import treecorr
import fitsio as fi


from mbii.readsubhalo import *
from mbii.properties import *


groups = mbii.mbdb.groups()

data1=fi.FITS('/home/ssamurof/massive_black_ii/halo_positions.fits')[-1].read()
reload(mbii.mbdb)
groups=mbii.mbdb.groups()
r, w = groups.calc_gg_all(data1,data1, save=True)
