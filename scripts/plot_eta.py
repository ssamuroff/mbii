import fitsio as fi
import treecorr
import numpy as np 
import pylab as plt
plt.switch_backend('agg')

baryons1=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat.fits')['baryons'][:]
baryons2=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat_reduced.fits')['baryons'][:]
dm1=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat.fits')['dm'][:]
dm2=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat_reduced.fits')['dm'][:]

mask1 = baryons1['c3']!=0 
cat1 = treecorr.Catalog(x=baryons1['x'][mask1], y=baryons1['y'][mask1], z=baryons1['z'][mask1], a=baryons1['c1'][mask1], b=baryons1['c2'][mask1], c=baryons1['c3'][mask1])
mask2 = baryons2['c3']!=0
cat2 = treecorr.Catalog(x=baryons2['x'][mask2], y=baryons2['y'][mask2], z=baryons2['z'][mask2], a=baryons2['c1'][mask2], b=baryons2['c2'][mask2], c=baryons2['c3'][mask2])

mask1 = dm1['c3']!=0 
cat3 = treecorr.Catalog(x=dm1['x'][mask1], y=dm1['y'][mask1], z=dm1['z'][mask1], a=dm1['c1'][mask1], b=dm1['c2'][mask1], c=dm1['c3'][mask1])
mask2 = dm2['c3']!=0
cat4 = treecorr.Catalog(x=dm2['x'][mask2], y=dm2['y'][mask2], z=dm2['z'][mask2], a=dm2['c1'][mask2], b=dm2['c2'][mask2], c=dm2['c3'][mask2])


dm_basic = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
dm_reduced = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
star_basic = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
star_reduced = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
star_basic.process(cat1,cat1)
star_reduced.process(cat2,cat2)
dm_basic.process(cat3,cat3)
dm_reduced.process(cat4,cat4)

plt.close() ; 
plt.plot(np.exp(dm_reduced.logr)/1e3, dm_reduced.xi,color='k', ls='--', label=r'DM, Reduced Inertia Tensor')
plt.plot(np.exp(dm_basic.logr)/1e3, dm_basic.xi,color='k', ls=':', label=r'DM, Basic Inertia Tensor')
plt.plot(np.exp(star_reduced.logr)/1e3, star_reduced.xi,color='purple', ls='--', label=r'Stars, Reduced Inertia Tensor')
plt.plot(np.exp(star_basic.logr)/1e3, star_basic.xi,color='purple', ls=':', label=r'Stars, Basic Inertia Tensor')
plt.axhline(0,color='k')
plt.legend(loc='upper right')

plt.ylim(1e-3,1)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('$\eta_e (r)$', fontsize=18)
plt.xlabel('Comoving Separation $r$ / $h^{-1}$ Mpc', fontsize=18)
plt.subplots_adjust(left=0.16,top=0.98,bottom=0.16, right=0.98)
plt.savefig('/home/ssamurof/etae_redit_vs_it-log-v0.png')
