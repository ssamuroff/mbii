import numpy as np
import yaml
import argparse
import mbii.mbdb as mb
import sys 
sys.path.append('/home/ssamurof/.local/lib/python2.7/site-packages/')
sys.path.append('/home/ssamurof/.local/lib/python2.6/site-packages/')

config = yaml.load(open(sys.argv[1]))
grp = mb.halos(fatal_errors=config['fatal_errors'])

dm = config['dm_shapes']
st = config['star_shapes']
nthresh = config['occupation_threshold']

data = grp.compile_data(dm_shapes=dm, star_shapes=st, nmin=nthresh)
