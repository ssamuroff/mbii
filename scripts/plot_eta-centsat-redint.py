import fitsio as fi
import treecorr
import numpy as np 
import pylab as plt
plt.switch_backend('agg')

#baryons2=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat-nthreshold5.fits')['baryons'][:]
baryons1=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat_reduced-nthreshold5.fits')['baryons'][:]
#dm2=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat-nthreshold5.fits')['dm'][:]
dm1=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat_reduced-nthreshold5.fits')['dm'][:]

select = (dm1['npart']>1000) & (baryons1['npart']>1000) & (np.isfinite(baryons1['x']) & np.isfinite(baryons1['y']) & np.isfinite(baryons1['z'])) & (baryons1['x']<100000) & (baryons1['y']<100000) & (baryons1['z']<100000)   
cflag=fi.FITS('/home/ssamurof/subhalo_central_flags.fits')[1].read()['central1'].astype(bool)

mask1 = (baryons1['c3']!=0) & select & cflag 
cat1 = treecorr.Catalog(x=baryons1['x'][mask1], y=baryons1['y'][mask1], z=baryons1['z'][mask1], a=baryons1['c1'][mask1], b=baryons1['c2'][mask1], c=baryons1['c3'][mask1])
mask2 = (baryons1['c3']!=0) & select & np.invert(cflag)
cat2 = treecorr.Catalog(x=baryons1['x'][mask2], y=baryons1['y'][mask2], z=baryons1['z'][mask2], a=baryons1['c1'][mask2], b=baryons1['c2'][mask2], c=baryons1['c3'][mask2])

mask1 = (dm1['c3']!=0) & select & cflag
cat3 = treecorr.Catalog(x=dm1['x'][mask1], y=dm1['y'][mask1], z=dm1['z'][mask1], a=dm1['c1'][mask1], b=dm1['c2'][mask1], c=dm1['c3'][mask1])
mask2 = (dm1['c3']!=0) & select & np.invert(cflag)
cat4 = treecorr.Catalog(x=dm1['x'][mask2], y=dm1['y'][mask2], z=dm1['z'][mask2], a=dm1['c1'][mask2], b=dm1['c2'][mask2], c=dm1['c3'][mask2])


dm_ss = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
dm_cc = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
dm_cs = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
star_ss = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
star_cc = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
star_cs = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
star_cc.process(cat1,cat1)
star_cs.process(cat1,cat2)
star_ss.process(cat2,cat2)
dm_ss.process(cat3,cat3)
dm_cs.process(cat3,cat4)
dm_cc.process(cat4,cat4)

plt.close() ; 
sel=(dm_cc.weight>=100)
plt.plot(np.exp(dm_cc.logr[sel])/1e3, dm_cc.xi[sel],color='purple', ls='--', label=r'Matter, $cc$')
sel=(dm_ss.weight>=100)
plt.plot(np.exp(dm_ss.logr[sel])/1e3, dm_ss.xi[sel],color='purple', ls='-', label=r'Matter, $ss$')
sel=(dm_cs.weight>=100)
plt.plot(np.exp(dm_cs.logr[sel])/1e3, dm_cs.xi[sel],color='purple', ls=':', label=r'Matter, $sc$')
sel=(star_cc.weight>=100)
plt.plot(np.exp(star_cc.logr[sel])/1e3, star_cc.xi[sel],color='k', ls='--', label=r'Baryons, $cc$')
sel=(star_ss.weight>=100)
plt.plot(np.exp(star_ss.logr[sel])/1e3, star_ss.xi[sel],color='k', ls='-', label=r'Baryons, $ss$')
sel=(star_cs.weight>=100)
plt.plot(np.exp(star_cs.logr[sel])/1e3, star_cs.xi[sel],color='k', ls=':', label=r'Baryons, $sc$')
plt.axhline(0,color='k')
plt.legend(loc='upper right')

plt.ylim(1e-3,1)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('$\omega_e (r)$', fontsize=18)
plt.xlabel('Comoving Separation $r$ / $h^{-1}$ Mpc', fontsize=18)
plt.subplots_adjust(left=0.16,top=0.98,bottom=0.16, right=0.98)
plt.savefig('/home/ssamurof/etae_redit_vs_reducedit-log-satcent.png')
