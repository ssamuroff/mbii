import fitsio as fi
import treecorr
import numpy as np 
import pylab as plt
plt.switch_backend('agg')

path = '/physics2/ssamurof/massive_black_ii/cats/base_subhalo_shapes-v2.fits'
baryons1=fi.FITS(path)[1].read()


#select = (dm1['npart']>1000) & (baryons1['npart']>1) & (np.isfinite(baryons1['x']) & np.isfinite(baryons1['y']) & np.isfinite(baryons1['z'])) & (baryons1['x']<100000) & (baryons1['y']<100000) & (baryons1['z']<100000) & (baryons1['x']>0) & (baryons1['y']>0) & (baryons1['z']>0)
cflag = baryons1['central']

mask1 = (cflag==1)
cat1 = treecorr.Catalog(x=baryons1['x'][mask1], y=baryons1['y'][mask1], z=baryons1['z'][mask1], a=baryons1['a1'][mask1], b=baryons1['a2'][mask1], c=baryons1['a3'][mask1])
mask2 = (cflag!=1)
cat2 = treecorr.Catalog(x=baryons1['x'][mask2], y=baryons1['y'][mask2], z=baryons1['z'][mask2], a=baryons1['a1'][mask2], b=baryons1['a2'][mask2], c=baryons1['a3'][mask2])

mask1 = (cflag==1)
cat3 = treecorr.Catalog(x=baryons1['x'][mask1], y=baryons1['y'][mask1], z=baryons1['z'][mask1], a=baryons1['a1_dm'][mask1], b=baryons1['a2_dm'][mask1], c=baryons1['a3_dm'][mask1])
mask2 = (cflag!=1)
cat4 = treecorr.Catalog(x=baryons1['x'][mask2], y=baryons1['y'][mask2], z=baryons1['z'][mask2], a=baryons1['a1_dm'][mask2], b=baryons1['a2_dm'][mask2], c=baryons1['a3_dm'][mask2])


dm_ss = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
dm_cc = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
dm_cs = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
star_ss = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
star_cc = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
star_cs = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
star_cc.process(cat1,cat1)
star_cs.process(cat1,cat2)
star_ss.process(cat2,cat2)
dm_cc.process(cat3,cat3)
dm_cs.process(cat3,cat4)
dm_ss.process(cat4,cat4)

plt.close() ; 
sel=(dm_cc.weight>=100)
plt.plot(np.exp(dm_cc.logr[sel]), dm_cc.xi[sel],color='purple', ls='--', label=r'Matter, $cc$')
sel=(dm_ss.weight>=100)
plt.plot(np.exp(dm_ss.logr[sel]), dm_ss.xi[sel],color='purple', ls='-', label=r'Matter, $ss$')
sel=(dm_cs.weight>=100)
plt.plot(np.exp(dm_cs.logr[sel]), dm_cs.xi[sel],color='purple', ls=':', label=r'Matter, $sc$')
sel=(star_cc.weight>=100)
plt.plot(np.exp(star_cc.logr[sel]), star_cc[sel].xi,color='k', ls='--', label=r'Baryons, $cc$')
sel=(star_ss.weight>=100)
plt.plot(np.exp(star_ss.logr[sel]), star_ss.xi[sel],color='k', ls='-', label=r'Baryons, $ss$')
sel=(star_cs.weight>=100)
plt.plot(np.exp(star_cs.logr[sel]), star_cs.xi[sel],color='k', ls=':', label=r'Baryons, $sc$')
plt.axhline(0,color='k')
plt.legend(loc='upper right')

plt.ylim(1e-3,1)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('$\omega (r)$', fontsize=18)
plt.xlabel('Comoving Separation $r$ / $h^{-1}$ Mpc', fontsize=18)
plt.subplots_adjust(left=0.16,top=0.98,bottom=0.16, right=0.98)
plt.savefig('/physics2/ssamurof/massive_black_ii/plots/2pt/omega_base-csdecomp.png')
