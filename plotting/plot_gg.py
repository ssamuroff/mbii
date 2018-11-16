import numpy as np
import pylab as plt
import scipy.optimize as opt
plt.style.use('y1a1')
plt.switch_backend('agg')

central_name='spatial_central'
gg = {}
gg_sym = {}
base = '/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/'
gg[(0,0)] = np.loadtxt(base+'fiducial/20bins/gg_corr_%d%d.txt'%(0,0)).T
gg_sym[(0,0)] = np.loadtxt(base+'symmetrised/20bins/gg_corr_%d%d.txt'%(0,0)).T

delta = 10**-0.06
fig = plt.figure()
ax = fig.add_subplot(2, 1, 1)
plt.errorbar(gg[0,0][0][:-2], gg[0,0][1][:-2], yerr=gg[0,0][-1][:-2], marker='.', linestyle='none', capsize=0.2, color='purple', label='Unsymmetrised')
plt.errorbar(gg_sym[0,0][0][:-2]/delta, abs(gg_sym[0,0][1][:-2]), yerr=gg_sym[0,0][-1][:-2], marker='D', linestyle='none', capsize=0.2, color='k', markersize=3.5, label='Symmetrised')
plt.xscale('log')
plt.yscale('log')
ax.set_yscale('log', nonposy='clip')
leg=plt.legend(loc='upper right', fontsize=15)
plt.ylabel(r'$\xi_{gg}$', fontsize=18)
plt.yticks([1e-1,1e0,10], ['$0.1$', '$1.0$', '$10$'])
leg.get_frame().set_alpha(0.75)
plt.xlim(0.1, 33.)
plt.ylim(0.02, 160)

plt.xticks([0.1,1,10],visible=False)
plt.yticks(fontsize=16)

ax = fig.add_subplot(2, 1, 2)
plt.axhline(0, color='k', ls=':')
plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=18)
plt.ylabel('Frac. Residual $f_{gg}$', fontsize=18)
plt.ylim(-0.1,0.21)
plt.xscale('log')
plt.xlim(0.1, 33.)
plt.xticks([0.1,1,10],['$0.1$','$1.0$','$10$'],fontsize=16)
plt.yticks([-0.1,0,0.1,0.2], fontsize=16)
fgg = (gg[0,0][1]-gg_sym[0,0][1])/gg[0,0][1]
plt.plot(gg[0,0][0][:-2], fgg[:-2], linestyle='--', lw=2., color='purple')


plt.subplots_adjust(hspace=0, left=0.155, bottom=0.148,top=0.98)
plt.savefig('/Users/hattifattener/Documents/ias/mbii/draft/figures/fracres-panels-gg-mbii-symmetrised.pdf')
plt.close()
