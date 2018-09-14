import numpy as np
import scipy.optimize as opt
import pylab as plt
plt.switch_backend('agg')
plt.style.use('y1a1')

snapshot = 85

def func(x, a, b):
	return 10**b * x**a


gg_massive_black = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/fiducial/20bins/gg_corr_00.txt').T
gg_massive_black_sym = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/symmetrised/20bins/gg_corr_00.txt').T
gg_illustris = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/illustris/gg_corr_00.txt').T
gg_illustris_sym = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/illustris/symmetrised/gg_corr_00.txt').T

fig=plt.figure()
ax=fig.add_subplot(211)
delta = 10**-0.06

plt.errorbar(gg_massive_black[0]*delta, gg_massive_black[1], yerr=gg_massive_black[-1], color='purple', marker='D', mec='purple', linestyle='none', capsize=0.2, markersize=3.5, label='Unsymmetrised (MassiveBlack-II)')
plt.errorbar(gg_massive_black_sym[0], gg_massive_black_sym[1], yerr=gg_massive_black_sym[-1], color='royalblue', mec='royalblue', marker='D', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Symmetrised (MassiveBlack-II)')
plt.errorbar(gg_illustris[0]*delta, gg_illustris[1], yerr=gg_illustris[-1], color='plum', marker='^', mec='plum', linestyle='none', capsize=0.2, markersize=3.5, label='Unsymmetrised (Illustris)')
plt.errorbar(gg_illustris_sym[0], gg_illustris_sym[1], yerr=gg_illustris_sym[-1], color='pink', mec='pink', marker='^', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Symmetrised (Illustris)')
ax.set_yscale('log', nonposy="clip")
plt.legend(loc='lower left', fontsize=12)
plt.xlim(0.07,33.)
plt.ylabel(r'$\xi_{gg}$', fontsize=14)
#plt.ylim(1e-5,1e-1)
#plt.yticks([1e-4, 1e-3,1e-2,1e-1], fontsize=13)
plt.xticks(visible=False)
plt.xscale('log')

ax=fig.add_subplot(212)
delta = 10**-0.06

p_gg,c_gg=opt.curve_fit(func, gg_massive_black_sym[0], gg_massive_black_sym[1], sigma=gg_massive_black_sym[-1])
gg_smooth_mb = func(gg_massive_black[0], p_gg[0], p_gg[1])
p_gg,c_gg=opt.curve_fit(func, gg_illustris_sym[0], gg_illustris_sym[1], sigma=gg_illustris_sym[-1])
gg_smooth_il = func(gg_illustris[0], p_gg[0], p_gg[1])

fm = ( gg_massive_black[1] -  gg_massive_black_sym[1])/ gg_massive_black[1]
fi = ( gg_illustris[1] -  gg_illustris_sym[1])/ gg_illustris[1]

plt.plot(gg_massive_black[0][:-2], fm[:-2], color='purple', linestyle='-', lw=2.0, label='MassiveBlack-II')
plt.plot(gg_illustris[0], fi, color='pink', linestyle='--', lw=2.0, label='Illustris')
plt.legend(loc='upper right', fontsize=12)
plt.xlim(0.07,33.)
plt.ylim(-0.1,0.2)
#plt.yticks([1e-4, 1e-3,1e-2,1e-1], fontsize=13)
plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=14)
plt.ylabel(r'Frac. Residual $f_{gg}$', fontsize=14)
plt.xscale('log')
plt.axhline(0, color='k', ls=':')

plt.subplots_adjust(hspace=0, wspace=0, left=0.14, right=0.95, top=0.98, bottom=0.14)
plt.savefig('/Users/hattifattener/Documents/ias/mbii/draft/figures/gg-symmetrised-split-simcompare.pdf')




plt.close()



ed_massive_black = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/fiducial/ed_corr_00.txt').T
ed_massive_black_sym = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/symmetrised/ed_corr_00.txt').T
ed_illustris = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/illustris/ed_corr_00_splitbyspatial_central.txt').T
ed_illustris_sym = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/illustris/symmetrised/ed_corr_00_splitbyspatial_central.txt').T

fig=plt.figure()
ax=fig.add_subplot(211)
delta = 10**-0.06

plt.errorbar(ed_massive_black[0]*delta, ed_massive_black[1], yerr=ed_massive_black[-1], color='purple', marker='D', mec='purple', linestyle='none', capsize=0.2, markersize=3.5, label='Unsymmetrised (MassiveBlack-II)')
plt.errorbar(ed_massive_black_sym[0], ed_massive_black_sym[1], yerr=ed_massive_black_sym[-1], color='royalblue', mec='royalblue', marker='D', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Symmetrised (MassiveBlack-II)')
plt.errorbar(ed_illustris[0]*delta, ed_illustris[1], yerr=ed_illustris[-1], color='plum', marker='^', mec='plum', linestyle='none', capsize=0.2, markersize=3.5, label='Unsymmetrised (Illustris)')
plt.errorbar(ed_illustris_sym[0], ed_illustris_sym[1], yerr=ed_illustris_sym[-1], color='pink', mec='pink', marker='^', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Symmetrised (Illustris)')
ax.set_yscale('log', nonposy="clip")
plt.legend(loc='lower left', fontsize=12)
plt.xlim(0.07,33.)
plt.ylabel(r'ED', fontsize=14)
#plt.ylim(1e-5,1e-1)
#plt.yticks([1e-4, 1e-3,1e-2,1e-1], fontsize=13)
plt.xticks(visible=False)
plt.xscale('log')

ax=fig.add_subplot(212)
delta = 10**-0.06

p_ed,c_ed=opt.curve_fit(func, ed_massive_black_sym[0], ed_massive_black_sym[1], sigma=ed_massive_black_sym[2])
ed_smooth_mb = func(ed_massive_black[0], p_ed[0], p_ed[1])
p_ed,c_ed=opt.curve_fit(func, ed_illustris_sym[0], ed_illustris_sym[1], sigma=ed_illustris_sym[2])
ed_smooth_il = func(ed_illustris[0], p_ed[0], p_ed[1])

fm = ( ed_massive_black[1] -  ed_massive_black_sym[1])/ ed_massive_black[1]
fi = ( ed_illustris[1] -  ed_illustris_sym[1])/ ed_illustris[1]

plt.plot(ed_massive_black[0], fm, color='purple', linestyle='-', lw=2.0, label='MassiveBlack-II')
plt.plot(ed_illustris[0], fi, color='pink', linestyle='--', lw=2.0, label='Illustris')
plt.legend(loc='upper right', fontsize=12)
plt.xlim(0.07,33.)
#plt.ylim(1e-5,1e-1)
#plt.yticks([1e-4, 1e-3,1e-2,1e-1], fontsize=13)
plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=14)
plt.ylabel(r'Frac. Residual $f_\mathrm{ED}$', fontsize=14)
plt.xscale('log')
plt.axhline(0, color='k', ls=':')

plt.subplots_adjust(hspace=0, wspace=0, left=0.14, right=0.95, top=0.98, bottom=0.14)
plt.savefig('/Users/hattifattener/Documents/ias/mbii/draft/figures/ed-symmetrised-split-simcompare.pdf')






plt.close()

ee_massive_black = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/fiducial/ee_corr_00.txt').T
ee_massive_black_sym = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/symmetrised/ee_corr_00.txt').T
ee_illustris = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/illustris/ee_corr_00_splitbyspatial_central.txt').T
ee_illustris_sym = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/illustris/symmetrised/ee_corr_00_splitbyspatial_central.txt').T

fig=plt.figure()
ax=fig.add_subplot(211)
delta = 10**-0.06

plt.errorbar(ee_massive_black[0]*delta, ee_massive_black[1], yerr=ee_massive_black[-1], color='purple', marker='D', mec='purple', linestyle='none', capsize=0.2, markersize=3.5, label='Unsymmetrised (MassiveBlack-II)')
plt.errorbar(ee_massive_black_sym[0], ee_massive_black_sym[1], yerr=ee_massive_black_sym[-1], color='royalblue', mec='royalblue', marker='D', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Symmetrised (MassiveBlack-II)')
plt.errorbar(ee_illustris[0]*delta, ee_illustris[1], yerr=ee_illustris[-1], color='plum', marker='^', mec='plum', linestyle='none', capsize=0.2, markersize=3.5, label='Unsymmetrised (Illustris)')
plt.errorbar(ee_illustris_sym[0], ee_illustris_sym[1], yerr=ee_illustris_sym[-1], color='pink', mec='pink', marker='^', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Symmetrised (Illustris)')
ax.set_yscale('log', nonposy="clip")
#plt.legend(loc='lower left', fontsize=12)
plt.xlim(0.07,33.)
plt.ylabel(r'EE', fontsize=14)
#plt.ylim(1e-5,1e-1)
#plt.yticks([1e-4, 1e-3,1e-2,1e-1], fontsize=13)
plt.xticks(visible=False)
plt.xscale('log')

ax=fig.add_subplot(212)
delta = 10**-0.06

p_ee,c_ee=opt.curve_fit(func, ee_massive_black_sym[0][ee_massive_black_sym[1]>0], ee_massive_black_sym[1][ee_massive_black_sym[1]>0], sigma=ee_massive_black_sym[2][ee_massive_black_sym[1]>0])
ee_smooth_mb = func(ee_massive_black[0], p_ee[0], p_ee[1])
p_ee,c_ee=opt.curve_fit(func, ee_illustris_sym[0][ee_illustris_sym[1]>0], ee_illustris_sym[1][ee_illustris_sym[1]>0], sigma=ee_illustris_sym[2][ee_illustris_sym[1]>0])
ee_smooth_il = func(ee_illustris[0], p_ee[0], p_ee[1])

fm = ( ee_massive_black[1] -  ee_massive_black_sym[1])/ ee_smooth_mb
fi = ( ee_illustris[1] -  ee_illustris_sym[1])/ ee_smooth_il

plt.plot(ee_massive_black[0], fm, color='purple', linestyle='-', lw=2.0, label='MassiveBlack-II')
plt.plot(ee_illustris[0], fi, color='pink', linestyle='--', lw=2.0, label='Illustris')
#plt.legend(loc='upper right', fontsize=12)
plt.xlim(0.07,33.)
plt.ylim(-0.5,2.5)
#plt.yticks([1e-4, 1e-3,1e-2,1e-1], fontsize=13)
plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=14)
plt.ylabel(r'Frac. Residual $f_\mathrm{EE}$', fontsize=14)
plt.xscale('log')
plt.axhline(0, color='k', ls=':')

plt.subplots_adjust(hspace=0, wspace=0, left=0.14, right=0.95, top=0.98, bottom=0.14)
plt.savefig('/Users/hattifattener/Documents/ias/mbii/draft/figures/ee-symmetrised-split-simcompare.pdf')




plt.close()

gi_massive_black = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/fiducial/GIplus_proj_corr_00.txt').T
gi_massive_black_sym = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/symmetrised/GIplus_proj_corr_00.txt').T
gi_illustris = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/illustris/GIplus_proj_corr_00_splitbyspatial_central.txt').T
gi_illustris_sym = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/illustris/symmetrised/GIplus_proj_corr_00_splitbyspatial_central.txt').T

fig=plt.figure()
ax=fig.add_subplot(211)
delta = 10**-0.06

plt.errorbar(gi_massive_black[0]*delta, gi_massive_black[1], yerr=gi_massive_black[-1], color='purple', marker='D', mec='purple', linestyle='none', capsize=0.2, markersize=3.5, label='Unsymmetrised (MassiveBlack-II)')
plt.errorbar(gi_massive_black_sym[0], gi_massive_black_sym[1], yerr=gi_massive_black_sym[-1], color='royalblue', mec='royalblue', marker='D', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Symmetrised (MassiveBlack-II)')
plt.errorbar(gi_illustris[0]*delta, gi_illustris[1], yerr=gi_illustris[-1], color='plum', marker='^', mec='plum', linestyle='none', capsize=0.2, markersize=3.5, label='Unsymmetrised (Illustris)')
plt.errorbar(gi_illustris_sym[0], gi_illustris_sym[1], yerr=gi_illustris_sym[-1], color='pink', mec='pink', marker='^', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Symmetrised (Illustris)')
ax.set_yscale('log', nonposy="clip")
#plt.legend(loc='lower left', fontsize=12)
plt.xlim(0.07,33.)
plt.ylabel(r'$w_{g+}$', fontsize=14)
#plt.ylim(1e-5,1e-1)
#plt.yticks([1e-4, 1e-3,1e-2,1e-1], fontsize=13)
plt.xticks(visible=False)
plt.xscale('log')

ax=fig.add_subplot(212)
delta = 10**-0.06

p_gi,c_gi=opt.curve_fit(func, gi_massive_black_sym[0], gi_massive_black_sym[1], sigma=gi_massive_black_sym[2])
gi_smooth_mb = func(gi_massive_black[0], p_gi[0], p_gi[1])
p_gi,c_gi=opt.curve_fit(func, gi_illustris_sym[0], gi_illustris_sym[1], sigma=gi_illustris_sym[2])
gi_smooth_il = func(gi_illustris[0], p_gi[0], p_gi[1])

fm = ( gi_massive_black[1] -  gi_massive_black_sym[1])/ gi_massive_black[1]
fi = ( gi_illustris[1] -  gi_illustris_sym[1])/ gi_illustris[1]

plt.plot(gi_massive_black[0], fm, color='purple', linestyle='-', lw=2.0, label='MassiveBlack-II')
plt.plot(gi_illustris[0], fi, color='pink', linestyle='--', lw=2.0, label='Illustris')
#plt.legend(loc='upper right', fontsize=12)
plt.xlim(0.07,33.)
#plt.ylim(-0.5,2.5)
#plt.yticks([1e-4, 1e-3,1e-2,1e-1], fontsize=13)
plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=14)
plt.ylabel(r'Frac. Residual $f_{g+}$', fontsize=14)
plt.xscale('log')
plt.axhline(0, color='k', ls=':')

plt.subplots_adjust(hspace=0, wspace=0, left=0.14, right=0.95, top=0.98, bottom=0.14)
plt.savefig('/Users/hattifattener/Documents/ias/mbii/draft/figures/wgp-symmetrised-split-simcompare.pdf')



plt.close()

ii_massive_black = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/fiducial/IIplus_proj_corr_00.txt').T
ii_massive_black_sym = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/symmetrised/IIplus_proj_corr_00.txt').T
ii_illustris = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/illustris/IIplus_proj_corr_00_splitbyspatial_central.txt').T
ii_illustris_sym = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/illustris/symmetrised/IIplus_proj_corr_00_splitbyspatial_central.txt').T

fig=plt.figure()
ax=fig.add_subplot(211)
delta = 10**-0.06

plt.errorbar(ii_massive_black[0]*delta, ii_massive_black[1], yerr=ii_massive_black[-1], color='purple', marker='D', mec='purple', linestyle='none', capsize=0.2, markersize=3.5, label='Unsymmetrised (MassiveBlack-II)')
plt.errorbar(ii_massive_black_sym[0], ii_massive_black_sym[1], yerr=ii_massive_black_sym[-1], color='royalblue', mec='royalblue', marker='D', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Symmetrised (MassiveBlack-II)')
plt.errorbar(ii_illustris[0]*delta, ii_illustris[1], yerr=ii_illustris[-1], color='plum', marker='^', mec='plum', linestyle='none', capsize=0.2, markersize=3.5, label='Unsymmetrised (Illustris)')
plt.errorbar(ii_illustris_sym[0], ii_illustris_sym[1], yerr=ii_illustris_sym[-1], color='pink', mec='pink', marker='^', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Symmetrised (Illustris)')
ax.set_yscale('log', nonposy="clip")
#plt.legend(loc='lower left', fontsize=12)
plt.xlim(0.07,33.)
plt.ylabel(r'$w_{++}$', fontsize=14)
#plt.ylim(1e-5,1e-1)
#plt.yticks([1e-4, 1e-3,1e-2,1e-1], fontsize=13)
plt.xticks(visible=False)
plt.xscale('log')

ax=fig.add_subplot(212)
delta = 10**-0.06

p_ii,c_ii=opt.curve_fit(func, ii_massive_black_sym[0], ii_massive_black_sym[1], sigma=ii_massive_black_sym[2])
ii_smooth_mb = func(ii_massive_black[0], p_ii[0], p_ii[1])
p_ii,c_ii=opt.curve_fit(func, ii_illustris_sym[0], ii_illustris_sym[1], sigma=ii_illustris_sym[2])
ii_smooth_il = func(gi_illustris[0], p_ii[0], p_ii[1])

fm = ( ii_massive_black[1] -  ii_massive_black_sym[1])/ ii_massive_black[1]
fi = ( ii_illustris[1] -  ii_illustris_sym[1])/ ii_illustris[1]

plt.plot(ii_massive_black[0], fm, color='purple', linestyle='-', lw=2.0, label='MassiveBlack-II')
plt.plot(ii_illustris[0], fi, color='pink', linestyle='--', lw=2.0, label='Illustris')
#plt.legend(loc='upper right', fontsize=12)
plt.xlim(0.07,33.)
#plt.ylim(-0.5,2.5)
#plt.yticks([1e-4, 1e-3,1e-2,1e-1], fontsize=13)
plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=14)
plt.ylabel(r'Frac. Residual $f_{++}$', fontsize=14)
plt.xscale('log')
plt.axhline(0, color='k', ls=':')

plt.subplots_adjust(hspace=0, wspace=0, left=0.14, right=0.95, top=0.98, bottom=0.14)
plt.savefig('/Users/hattifattener/Documents/ias/mbii/draft/figures/wpp-symmetrised-split-simcompare.pdf')




gg_massive_black = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/fiducial/illustris_binning/gg_corr_00_splitbyspatial_central.txt').T
gg_massive_black_sym = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/symmetrised/illustris_binning/gg_corr_00_splitbyspatial_central.txt').T
ed_massive_black = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/fiducial/illustris_binning/ed_corr_00_splitbyspatial_central.txt').T
ed_massive_black_sym = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/symmetrised/illustris_binning/ed_corr_00_splitbyspatial_central.txt').T
ee_massive_black = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/fiducial/illustris_binning/ee_corr_00_splitbyspatial_central.txt').T
ee_massive_black_sym = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/symmetrised/illustris_binning/ee_corr_00_splitbyspatial_central.txt').T
gi_massive_black = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/fiducial/illustris_binning/GIplus_proj_corr_00_splitbyspatial_central.txt').T
gi_massive_black_sym = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/symmetrised/illustris_binning/GIplus_proj_corr_00_splitbyspatial_central.txt').T
ii_massive_black = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/fiducial/illustris_binning/IIplus_proj_corr_00_splitbyspatial_central.txt').T
ii_massive_black_sym = np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/symmetrised/illustris_binning/IIplus_proj_corr_00_splitbyspatial_central.txt').T


plt.close()

fig=plt.figure()
ax=fig.add_subplot(111)
delta = 10**-0.06

p_ee,c_ee=opt.curve_fit(func, ee_massive_black_sym[0][ee_massive_black_sym[1]>0], ee_massive_black_sym[1][ee_massive_black_sym[1]>0], sigma=ee_massive_black_sym[2][ee_massive_black_sym[1]>0])
ee_massive_black_smooth = func(ee_massive_black[0], p_ee[0], p_ee[1])
p_ee,c_ee=opt.curve_fit(func, ee_illustris_sym[0][ee_illustris_sym[1]>0], ee_illustris_sym[1][ee_illustris_sym[1]>0], sigma=ee_illustris_sym[2][ee_illustris_sym[1]>0])
ee_illustris_smooth = func(ee_illustris[0], p_ee[0], p_ee[1])

p_ed,c_ed=opt.curve_fit(func, ed_massive_black_sym[0][ed_massive_black_sym[1]>0], ed_massive_black_sym[1][ed_massive_black_sym[1]>0], sigma=ed_massive_black_sym[2][ed_massive_black_sym[1]>0])
ed_massive_black_smooth = func(ed_massive_black[0], p_ed[0], p_ed[1])
p_ed,c_ed=opt.curve_fit(func, ed_illustris_sym[0][ed_illustris_sym[1]>0], ed_illustris_sym[1][ed_illustris_sym[1]>0], sigma=ed_illustris_sym[2][ed_illustris_sym[1]>0])
ed_illustris_smooth = func(ed_illustris[0], p_ed[0], p_ed[1])


p_gg,c_gg=opt.curve_fit(func, gg_massive_black_sym[0][gg_massive_black_sym[1]>0], gg_massive_black_sym[1][gg_massive_black_sym[1]>0], sigma=gg_massive_black_sym[2][gg_massive_black_sym[1]>0])
gg_massive_black_smooth = func(gg_massive_black[0], p_gg[0], p_gg[1])
p_gg,c_gg=opt.curve_fit(func, gg_illustris_sym[0][gg_illustris_sym[1]>0], gg_illustris_sym[1][gg_illustris_sym[1]>0], sigma=gg_illustris_sym[2][gg_illustris_sym[1]>0])
gg_illustris_smooth = func(gg_illustris[0], p_gg[0], p_gg[1])


p_gi,c_gi=opt.curve_fit(func, gi_massive_black_sym[0][gi_massive_black_sym[1]>0], gi_massive_black_sym[1][gi_massive_black_sym[1]>0], sigma=gi_massive_black_sym[2][gi_massive_black_sym[1]>0])
gi_massive_black_smooth = func(gi_massive_black[0], p_gi[0], p_gi[1])
p_gi,c_gi=opt.curve_fit(func, gi_illustris_sym[0][gi_illustris_sym[1]>0], gi_illustris_sym[1][gi_illustris_sym[1]>0], sigma=gi_illustris_sym[2][gi_illustris_sym[1]>0])
gi_illustris_smooth = func(gi_illustris[0], p_gi[0], p_gi[1])


p_ii,c_ii=opt.curve_fit(func, ii_massive_black_sym[0][ii_massive_black_sym[1]>0], ii_massive_black_sym[1][ii_massive_black_sym[1]>0], sigma=ii_massive_black_sym[2][ii_massive_black_sym[1]>0])
ii_massive_black_smooth = func(ii_massive_black[0], p_ii[0], p_ii[1])
p_ii,c_gg=opt.curve_fit(func, ii_illustris_sym[0][ii_illustris_sym[1]>0], ii_illustris_sym[1][ii_illustris_sym[1]>0], sigma=ii_illustris_sym[2][ii_illustris_sym[1]>0])
ii_illustris_smooth = func(ii_illustris[0], p_ii[0], p_ii[1])





fgg = ( gg_massive_black[1] -  gg_massive_black_sym[1])/gg_massive_black_smooth - ( gg_illustris[1] -  gg_illustris_sym[1]) / gg_illustris_smooth
fee = ( ee_massive_black[1] -  ee_massive_black_sym[1])/ee_massive_black_smooth - ( ee_illustris[1] -  ee_illustris_sym[1]) / ee_illustris_smooth
fed = ( ed_massive_black[1] -  ed_massive_black_sym[1])/ed_massive_black_smooth - ( ed_illustris[1] -  ed_illustris_sym[1]) / ed_illustris_smooth
fgi = ( gi_massive_black[1] -  gi_massive_black_sym[1])/gi_massive_black_smooth - ( gi_illustris[1] -  gi_illustris_sym[1]) / gi_illustris_smooth
fii = ( ii_massive_black[1] -  ii_massive_black_sym[1])/ii_massive_black_smooth - ( ii_illustris[1] -  ii_illustris_sym[1]) / ii_illustris_smooth


plt.plot(gg_massive_black[0], fgg, color='k', linestyle=':', lw=2.0, label=r'$\delta_{g}\delta_{g}$')
plt.plot(ed_massive_black[0], fed, color='purple', linestyle='-', lw=2.0, label='ED')
plt.plot(ee_massive_black[0], fee, color='pink', linestyle='--', lw=2.0, label='EE')
plt.plot(gi_massive_black[0], fgi, color='plum', linestyle='-.', lw=2.0, label='$g+$')
plt.plot(gi_massive_black[0], fii, color='steelblue', linestyle=':', lw=2.0, label='$++$')
#plt.legend(loc='upper right', fontsize=12)
plt.xlim(0.07,33.)
#plt.ylim(-0.08,0.1)
#plt.yticks([1e-4, 1e-3,1e-2,1e-1], fontsize=13)
plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=14)
plt.ylabel(r'Change in Frac. Residual $\delta f_{\alpha}$', fontsize=16)
plt.xscale('log')
plt.axhline(0, color='k', ls=':')
plt.legend(loc='lower right')

plt.subplots_adjust(hspace=0, wspace=0, left=0.18, right=0.95, top=0.98, bottom=0.14)
plt.savefig('/Users/hattifattener/Documents/ias/mbii/draft/figures/all-symmetrised-split-simcompare.pdf')