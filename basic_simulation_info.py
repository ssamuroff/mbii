import astropy
import numpy as numpy

snapshots={85:0}
particle_mass = {'m':1.1e7, 'b':2.2e6}

cosmology = astropy.cosmology.FlatLambdaCDM(H0=70.2, Om0=0.275)

Lbox = 100 * 1e3 # In units of h^-1 kpc

cosmoparam = {'omega_m':0.275, 'sigma8':0.816, 'ns': 0.968, 'omega_b': 0.046, 'omega_de': 0.725, 'h': 0.701,'w': -1.0}