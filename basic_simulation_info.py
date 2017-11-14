import astropy
import numpy as numpy

snapshots={85:0}

cosmology = astropy.cosmology.FlatLambdaCDM(H0=70.2, Om0=0.275)

Lbox = 100 * 1e3 # In units of h^-1 kpc