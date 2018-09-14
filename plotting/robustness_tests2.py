import numpy as np
import pylab as plt
import scipy.optimize as opt
plt.style.use('y1a1')
plt.switch_backend('agg')

def func(x, a, b):
	return 10**b * x**a

def chi2(d, err):
	return sum(d*d/err/err)/len(d)

nsub=3
fac = np.sqrt(nsub**3)

snapshot=85
ed=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/draft/data/robustness/ed_corrvar_ref.txt').T
ed1=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/draft/data/robustness/ed_corrvar_altcent.txt').T
ed2=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/draft/data/robustness/ed_corrvar_galaxygalaxysym.txt').T
ed3=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/draft/data/robustness/ed_corrvar_dmshapes.txt').T

ee=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/draft/data/robustness/ee_corrvar_ref.txt').T
ee1=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/draft/data/robustness/ee_corrvar_altcent.txt').T
ee2=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/draft/data/robustness/ee_corrvar_galaxygalaxysym.txt').T
ee3=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/draft/data/robustness/ee_corrvar_dmshapes.txt').T

gi=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/draft/data/robustness/gi_plus_projected_corrvar_ref.txt').T
gi1=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/draft/data/robustness/gi_plus_projected_corrvar_altcent.txt').T
gi2=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/draft/data/robustness/gi_plus_projected_corrvar_galaxygalaxysym.txt').T
gi3=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/draft/data/robustness/gi_plus_projected_corrvar_dmshapes.txt').T

#ii=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/draft/data/robustness/ii_plus_projected_corrvar_ref.txt').T
ii1=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/draft/data/robustness/ii_plus_projected_corrvar_altcent.txt').T
ii2=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/draft/data/robustness/ii_plus_projected_corrvar_galaxygalaxysym.txt').T
ii3=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/draft/data/robustness/ii_plus_projected_corrvar_dmshapes.txt').T


fig=plt.figure()
delta = 10**-0.05


plt.close()
plt.subplot(411)
c1 = chi2(ed1[1],ed1[2]/fac)
c2 = chi2(ed2[1],ed2[2]/fac)
c3 = chi2(ed3[1],ed3[2]/fac)

plt.errorbar(ed1[0]*delta, ed1[1], yerr=ed1[2]/fac, linestyle='none', lw=1.0, ms=3.0, marker='D', mfc='forestgreen', color='forestgreen', mec='forestgreen', label=r'Alt. Central Dfn $\chi^2/\nu = %2.2f$'%c1)
plt.errorbar(ed2[0], ed2[1], yerr=ed2[2]/fac, linestyle='none', lw=1.0, marker='.', mfc='pink', color='pink', mec='pink', label=r'Galaxy-Galaxy Sym.  $\chi^2/\nu = %2.2f$'%c2)
plt.errorbar(ed3[0]/delta, ed3[1], yerr=ed3[2]/fac, linestyle='none', lw=1.0, ms=3.0, marker='D', mfc='none', color='royalblue', label=r'Subhalo Shapes  $\chi^2/\nu = %2.2f$'%c3)
#plt.fill_between(ed[0],-ed[1],ed[1],color='gray',alpha=0.1)
plt.xlim(0.1,33.)
plt.ylabel(r'$\delta f_\mathrm{ED}$', fontsize=16)
plt.xscale('log')
plt.axhline(0,color='k', ls=':')
plt.xticks(visible=False)
plt.legend(loc='lower left', fontsize=7)
plt.yticks([-0.4,-0.2,0.0], fontsize=12)
plt.ylim(-0.45,0.15)


plt.subplot(412)
c1 = chi2(ee1[1],ee1[2]/fac)
c2 = chi2(ee2[1],ee2[2]/fac)
c3 = chi2(ee3[1],ee3[2]/fac)

plt.errorbar(ee1[0]*delta, ee1[1], yerr=ee1[2]/fac, linestyle='none', lw=1.0, ms=3.0, marker='D', mfc='forestgreen', color='forestgreen', mec='forestgreen', label=r'$\chi^2/\nu = %2.2f$'%c1)
plt.errorbar(ee2[0], ee2[1], yerr=ee2[2]/fac, linestyle='none', lw=1.0, marker='.', mfc='pink', color='pink', mec='pink', label=r'$\chi^2/\nu = %2.2f$'%c2)
plt.errorbar(ee3[0]/delta, ee3[1], yerr=ee3[2]/fac, linestyle='none', lw=1.0, ms=3.0, marker='D', mfc='none', color='royalblue', label=r'$\chi^2/\nu = %2.2f$'%c3)
plt.xlim(0.1,33.)
plt.ylabel(r'$\delta f_\mathrm{EE}$', fontsize=16)
plt.xscale('log')
plt.axhline(0,color='k', ls=':')
plt.xticks(visible=False)
plt.legend(loc='upper left', fontsize=7)
plt.yticks([-1.2,-0.8,-0.4,0.0,0.,0.4,0.8], fontsize=12)
plt.ylim(-0.8,1.55)


plt.subplot(413)
c1 = chi2(gi1[1],gi1[2]/fac)
c2 = chi2(gi2[1],gi2[2]/fac)
c3 = chi2(gi3[1],gi3[2]/fac)

plt.errorbar(gi1[0]*delta, gi1[1], yerr=gi1[2]/fac, linestyle='none', lw=1.0, ms=3.0, marker='D', mfc='forestgreen', color='forestgreen', mec='forestgreen', label=r'$\chi^2/\nu = %2.2f$'%c1)
plt.errorbar(gi2[0], gi2[1], yerr=gi2[2]/fac, linestyle='none', lw=1.0, marker='.', mfc='pink', color='pink', mec='pink', label=r'$\chi^2/\nu = %2.2f$'%c2)
plt.errorbar(gi3[0]/delta, gi3[1], yerr=gi3[2]/fac, linestyle='none', lw=1.0, ms=3.0, marker='D', mfc='none', color='royalblue', label=r'$\chi^2/\nu = %2.2f$'%c3)
plt.xlim(0.1,33.)
plt.ylabel(r'$\delta f_\mathrm{g+}$', fontsize=16)
plt.xscale('log')
plt.axhline(0,color='k', ls=':')
plt.xticks(visible=False)
plt.legend(loc='upper left', fontsize=7)
plt.yticks([-1.2,-0.8,-0.4,0.0,0.,0.4,0.8], fontsize=12)
plt.ylim(-0.8,1.55)

plt.subplot(414)
c1 = chi2(ii1[1],ii1[2]/fac)
c2 = chi2(ii2[1],ii2[2]/fac)
c3 = chi2(ii3[1],ii3[2]/fac)

plt.errorbar(ii1[0]*delta, ii1[1], yerr=ii1[2]/fac, linestyle='none', lw=1.0, ms=3.0, marker='D', mfc='forestgreen', color='forestgreen', mec='forestgreen', label=r'$\chi^2/\nu = %2.2f$'%c1)
plt.errorbar(ii2[0], ii2[1], yerr=ii2[2]/fac, linestyle='none', lw=1.0, marker='.', mfc='pink', color='pink', mec='pink', label=r'$\chi^2/\nu = %2.2f$'%c2)
plt.errorbar(ii3[0]/delta, ii3[1], yerr=ii3[2]/fac, linestyle='none', lw=1.0, ms=3.0, marker='D', mfc='none', color='royalblue', label=r'$\chi^2/\nu = %2.2f$'%c3)
plt.xlim(0.1,33.)
plt.ylabel(r'$\delta f_\mathrm{++}$', fontsize=16)
plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=14)
plt.xscale('log')
plt.axhline(0,color='k', ls=':')
plt.xticks(visible=True)
plt.legend(loc='upper left', fontsize=7)
plt.yticks([-1.2,-0.8,-0.4,0.0,0.,0.4,0.8], fontsize=12)
plt.ylim(-0.8,1.55)


plt.subplots_adjust(hspace=0,wspace=0, bottom=0.145, left=0.15, right=0.65, top=0.97)
plt.savefig('/Users/hattifattener/Documents/ias/mbii/draft/figures/df-vs-r-robustness-v2.pdf')
