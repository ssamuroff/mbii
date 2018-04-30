import fitsio as fi
import treecorr
import numpy as np 
import mbii.lego_tools as util
import pylab as plt
plt.switch_backend('agg')

baryons1=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat-nthreshold5.fits')['baryons'][:]
baryons_symm = util.gather_mpi_output('/home/ssamurof/massive_black_ii/cats/symmetrised/*-masswtdmedian.fits', hdu='baryons')
dm1=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat-nthreshold5.fits')['dm'][:]

# The two catalogues hardcoded here should have identical selection functions
select = (dm1['npart']>1000) & (baryons1['npart']>1) & (np.isfinite(baryons1['x']) & np.isfinite(baryons1['y']) & np.isfinite(baryons1['z'])) & (baryons1['x']<100000) & (baryons1['y']<100000) & (baryons1['z']<100000) & (baryons1['x']>0) & (baryons1['y']>0) & (baryons1['z']>0)

mask1 = (baryons1['c3']!=0) & select 
cat1 = treecorr.Catalog(x=baryons1['x'][mask1], y=baryons1['y'][mask1], z=baryons1['z'][mask1], a=baryons1['c1'][mask1], b=baryons1['c2'][mask1], c=baryons1['c3'][mask1])
mask2 = (baryons_symm['c3']!=0) 
cat2 = treecorr.Catalog(x=baryons_symm['x'][mask2], y=baryons_symm['y'][mask2], z=baryons_symm['z'][mask2], a=baryons_symm['c1'][mask2], b=baryons_symm['c2'][mask2], c=baryons_symm['c3'][mask2])

star_asymm = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
star_symm = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
star_asymm.process(cat1,cat1, metric='Rperp')
star_symm.process(cat2,cat2, metric='Rperp')

plt.close() ; 
sel=(star_asymm.weight>=100)
plt.plot(np.exp(star_asymm.logr[sel])/1e3, star_asymm.xi[sel],color='purple', ls='-', label=r'Baryons')
sel=(star_symm.weight>=100)
plt.plot(np.exp(star_symm.logr[sel])/1e3, star_symm.xi[sel],color='k', ls='--', label=r'Baryons, symmetrised')

plt.axhline(0,color='k')
plt.legend(loc='upper right')

plt.ylim(0.5e-3,1)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('$\omega_e (r)$', fontsize=18)
plt.xlabel('Comoving Separation $r_p$ / $h^{-1}$ Mpc', fontsize=18)
plt.subplots_adjust(left=0.16,top=0.98,bottom=0.16, right=0.98)
plt.savefig('/home/ssamurof/etae_rperp_baseit-log-symm.png')
