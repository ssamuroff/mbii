import fitsio as fi
import treecorr
import numpy as np 
import pylab as plt
plt.switch_backend('agg')

baryons1=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat-nthreshold5.fits')['baryons'][:]
#baryons2=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat_reduced-nthreshold5.fits')['baryons'][:]
dm1=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat-nthreshold5.fits')['dm'][:]
#dm2=fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat_reduced-nthreshold5.fits')['dm'][:]

select = (dm1['npart']>1000) & (baryons1['npart']>1000) & (np.isfinite(baryons1['x']) & np.isfinite(baryons1['y']) & np.isfinite(baryons1['z'])) & (baryons1['x']<100000) & (baryons1['y']<100000) & (baryons1['z']<100000) & (baryons1['x']>0) & (baryons1['y']>0) & (baryons1['z']>0)
cflag=fi.FITS('/home/ssamurof/subhalo_central_flags.fits')[1].read()['central1'].astype(bool)

dr = np.sqrt((baryons1['x']-dm1['x'])**2 + (baryons1['y']-dm1['y'])**2 + (baryons1['z']-dm1['z'])**2)/1000



Mb = np.log10(baryons1['npart'] * 2.2e6)
M = np.log10(dm1['npart'] * 1.1e7)

plt.hist2d(M[select & cflag], dr[select & cflag], range=((9,12.5),(0,0.4)), bins=60)
plt.savefig('/home/ssamurof/massive_black_ii/plots/sanity/2dhist-central-bmoffset-mass.png')


plt.hist2d(Mb[select & cflag], dr[select & cflag], range=((7,12),(0,0.4)), bins=60)
plt.savefig('/home/ssamurof/massive_black_ii/plots/sanity/2dhist-central-bmoffset-barmass.png')

plt.hist2d(M[select & np.invert(cflag)], dr[select & np.invert(cflag)], range=((9,12.5),(0,0.4)), bins=60)
plt.savefig('/home/ssamurof/massive_black_ii/plots/sanity/2dhist-satellite-bmoffset-mass.png')


plt.hist2d(Mb[select & np.invert(cflag)], dr[select & np.invert(cflag)], range=((7,12),(0,0.4)), bins=60)
plt.savefig('/home/ssamurof/massive_black_ii/plots/sanity/2dhist-satellite-bmoffset-barmass.png')

exit()








Hc,xc = np.histogram(dr[select & cflag]/1000, range=(0,0.4), bins=65, normed=1)
Hs,xs = np.histogram(dr[select & np.invert(cflag)]/1000, range=(0,0.4), bins=65, normed=1)

x = (xs[:-1]+xs[1:])/2 

plt.fill_between(x, Hc, color='purple', alpha=0.2, lw=1.5, label='Centrals')
plt.plot(x, Hs, color='royalblue', lw=1.5, label='Satellites')
plt.legend(loc='upper right', fontsize=18)
plt.xlabel('Baryon - Matter offset $\delta r$ / Mpc $h^{-1}$', fontsize=18)
plt.yticks(visible=False)
plt.subplots_adjust(bottom=0.14)
plt.ylim(ymin=0)

plt.savefig('/home/ssamurof/massive_black_ii/plots/sanity/baryon-matter_offset-subhalos.png')
