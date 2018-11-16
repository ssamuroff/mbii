import numpy as np
from hankel import HankelTransform
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
import scipy.special as sps
import scipy.integrate as sint
import scipy.interpolate as spi
import pylab as plt

PgI={}
PII={}

k=np.logspace(-2,2,500)

colours = ['purple', 'pink', 'steelblue', 'plum']
redshift = {85: 0.06, 79: 0.3, 73 : 0.625, 68 : 1.0 }

for i, snapshot in enumerate([85, 79, 73, 68]):
  print i, snapshot
  # Load the measured correlation function
  gi=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/fiducial/GIplus_proj_corr_00.txt'%snapshot).T
  ii=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/fiducial/IIplus_proj_corr_00.txt'%snapshot).T

  # Initialise the transform
  ht4 = HankelTransform(nu=2, N=12000, h=0.005)
  ht4 = HankelTransform(nu=4, N=12000, h=0.005)
  spline = Spline(gi[0], gi[1], k=1)
  K = lambda x : spline(x)
  splineii = Spline(ii[0], ii[1], k=1)
  Kii = lambda x : splineii(x)

  # Evaluate the integral repeatedly at each of a set of k values 
  pk = np.array([ht4.transform(K, k0)[0] for k0 in k ])
  pkII = np.array([ 0.5 * (ht4.transform(Kii, k0)[0] + ht4.transform(Kii, k0)[0]) for k0 in k ])

  # Store the power spectrum for this redshift
  PgI[snapshot] = -1 * pk
  PII[snapshot] = pkII

  plt.plot(k, pk, color=colours[i], label='$z=%3.2f$'%redshift[snapshot])
  plt.plot(k, pkII, color=colours[i], ls='--')
  np.savetxt('pgI-massive_black_ii-snapshot%d.txt'%snapshot, [k,pk])
  np.savetxt('pII-massive_black_ii-snapshot%d.txt'%snapshot, [k,pkII])

plt.xlabel('Wavenumber $k$ / $h^{-1}$', fontsize=16)
plt.ylabel(r'$P_\mathrm{\delta_g I}(k)$', fontsize=16)
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-6,1e2)
plt.legend(loc='upper right')
plt.savefig('pgi-massive_black_ii.pdf')

