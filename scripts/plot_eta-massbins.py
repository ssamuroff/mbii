import fitsio as fi
import treecorr
import numpy as np 
import pylab as plt
plt.switch_backend('agg')


print 'Loading data'
baryons1=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat_reduced-nthreshold5.fits')['baryons'][:]
#baryons2=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat_reduced.fits')['baryons'][:]
dm1=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat_reduced-nthreshold5.fits')['dm'][:]
#dm2=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat_reduced.fits')['dm'][:]

select = (dm1['npart']>1000)
baryons1 = baryons1[select]
dm1 = dm1[select]

print 'Initialising catalogues'
Md = dm1['npart']*1.1e7

mask0 = (baryons1['c3']!=0) & (Md<1e11)
cat0 = treecorr.Catalog(x=baryons1['x'][mask0], y=baryons1['y'][mask0], z=baryons1['z'][mask0], a=baryons1['c1'][mask0], b=baryons1['c2'][mask0], c=baryons1['c3'][mask0])
mask1 = (baryons1['c3']!=0)
cat1 = treecorr.Catalog(x=baryons1['x'][mask1], y=baryons1['y'][mask1], z=baryons1['z'][mask1], a=baryons1['c1'][mask1], b=baryons1['c2'][mask1], c=baryons1['c3'][mask1])
mask2 = (baryons1['c3']!=0) & (Md>1e11)
cat2 = treecorr.Catalog(x=baryons1['x'][mask2], y=baryons1['y'][mask2], z=baryons1['z'][mask2], a=baryons1['c1'][mask2], b=baryons1['c2'][mask2], c=baryons1['c3'][mask2])

mask3 = (baryons1['c3']!=0 ) & (Md>1e12)
cat3 = treecorr.Catalog(x=baryons1['x'][mask3], y=baryons1['y'][mask3], z=baryons1['z'][mask3], a=baryons1['c1'][mask3], b=baryons1['c2'][mask3], c=baryons1['c3'][mask3])
mask4 = (baryons1['c3']!=0) & (Md>1e13)
cat4 = treecorr.Catalog(x=baryons1['x'][mask4], y=baryons1['y'][mask4], z=baryons1['z'][mask4], a=baryons1['c1'][mask4], b=baryons1['c2'][mask4], c=baryons1['c3'][mask4])



print 'Constructing correlations'
star0 = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
star1 = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
star2 = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
star3 = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
star4 = treecorr.NVCorrelation(min_sep=100, max_sep=1e5, nbins=16)
star0.process(cat0,cat0)
star1.process(cat1,cat1)
star2.process(cat2,cat2)
star3.process(cat3,cat3)
star4.process(cat4,cat4)

plt.close() ; 
plt.plot(np.exp(star1.logr)/1e3, star1.xi,color='purple', ls='-',lw=2, label=r'Baryons, all')
plt.plot(np.exp(star0.logr)/1e3, star0.xi,color='k', ls=':', label=r'Baryons, $M_\mathrm{m}<10^{11}$')
plt.plot(np.exp(star2.logr)/1e3, star2.xi,color='royalblue', ls=':', label=r'Baryons, $M_\mathrm{m}>10^{11}$')
plt.plot(np.exp(star3.logr)/1e3, star3.xi,color='royalblue', ls='--', label=r'Baryons, $M_\mathrm{m}>10^{12}$')
plt.plot(np.exp(star4.logr)/1e3, star4.xi,color='royalblue', ls='-.', label=r'Baryons, $M_\mathrm{m}>10^{13}$')
plt.axhline(0,color='k')
plt.legend(loc='lower left')

plt.ylim(1e-3,1)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('$\eta_e (r)$', fontsize=18)
plt.xlabel('Comoving Separation $r$ / $h^{-1}$ Mpc', fontsize=18)
plt.subplots_adjust(left=0.16,top=0.98,bottom=0.16, right=0.98)
plt.savefig('/home/ssamurof/etae_redit_vs_it-log-masscuts-v2.png')
