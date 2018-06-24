import argparse
import numpy as np
import yaml
import fitsio as fi
from halotools.mock_observables.alignments import ee_3d


pth='/Users/hattifattener/local/lib/python2.7/site-packages/mbii/config/fiducial_cat_hatti.yaml'
options = yaml.load(open(pth))
binning = options['2pt']['binning']
splitflag=options['2pt']['split']
data = fi.FITS(options['2pt']['shapes'])[-1].read()
name = options['2pt']['split']
print 'Dividing catalogue by %s'%name

mask = (data[name]>=options['2pt']['split_val'])

cat1 = data[mask]
cat2 = data[np.invert(mask)]
cat0 = data
data_sym = fi.FITS(options['2pt']['shapes'].replace('.fits', '-symmetrised.fits'))[-1].read()
cat1_sym = data_sym[mask]
cat2_sym = data_sym[np.invert(mask)]
cat0_sym = data_sym


avec2 = np.vstack((cat2_sym['a1'], cat2_sym['a2'], cat2_sym['a3'])).T
avec1 = np.vstack((cat2_sym['a1'], cat2_sym['a2'], cat2_sym['a3'])).T
pvec1 = np.vstack((cat2['x'], cat2['y'], cat2['z'])).T
pvec2 = np.vstack((cat2['x'], cat2['y'], cat2['z'])).T
rbins = np.logspace(np.log10(options['2pt']['rmin']), np.log10(options['2pt']['rmax']), options['2pt']['nbin'])

period=100.
ee = ee_3d(pvec2, avec2, pvec1, avec1, rbins, period=period, num_threads=1)
pvec1 = np.vstack((cat2_sym['x'], cat2_sym['y'], cat2_sym['z'])).T                 
pvec2 = np.vstack((cat2_sym['x'], cat2_sym['y'], cat2_sym['z'])).T
ee0 = ee_3d(pvec2, avec2, pvec1, avec1, rbins, period=period, num_threads=1)

dee=ee0-ee

plt.close()
plt.subplot(211)
plt.plot(x, ee0, color='purple', ls='--', label='Symmetrised', lw=2.5)
plt.plot(x, ee, color='plum', ls=':', label='Symmetrised Shapes', lw=2.5)
plt.xscale('log')
plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=18)
plt.ylim(1e-3,1e-1)
plt.ylabel('EE Correlation', fontsize=18)
plt.legend(loc='lower left')
plt.yscale('log')
plt.subplot(212)
plt.ylabel('Fractional Shift $\Delta$EE', fontsize=18)
plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=18)
plt.plot(x, dee/ee0, color='purple', ls='-',lw=2.5)
plt.xscale('log')
plt.subplots_adjust(hspace=0)
plt.xticks(fontsize=16)
#  [1e-1,1e0,1e1,1e2],['$0.1$','$1.0$', '$10.0$', '$100.0$'] ,fontsize=16)
plt.ylabel('Frac. Shift $\Delta$EE$/$EE', fontsize=18)
plt.yticks(fontsize=16)
 # [-0.01,0.01,0.03, 0.05],['$-0.01$','$0.01$', '$0.03$', '$0.05$'] ,fontsize=16)
plt.subplot(211)
plt.yticks([1e-6,1e-4,1e-2,1e0],['$10^{-6}$','$10^{-4}$', '$10^{-2}$', '$10^{0}$'] ,fontsize=16)
plt.xticks(visible=False)                                                            
plt.subplot(212)
plt.axhline(0,color='k', ls=':')
plt.savefig('/Users/hattifattener/Documents/ias/mbii/ee_corr-ss_binshift.pdf')