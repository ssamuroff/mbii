import numpy as np
import pylab as plt
import scipy.optimize as opt
plt.style.use('y1a1')
plt.switch_backend('agg')

def func(x, a, b):
	return 10**b * x**a


snapshot=85
ed=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/fiducial/ED_corr_00.txt'%snapshot).T
ed12=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/fiducial/12bins/ED_corr_00.txt'%snapshot).T
ed12_sym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/symmetrised/12bins/ED_corr_00.txt'%snapshot).T
ed_red=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/reduced/ED_corr_00_splitbyspatial_central.txt'%snapshot).T
ed_red_sym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/reduced/symmetrised/ED_corr_00_splitbyspatial_central.txt'%snapshot).T
ed_gsym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/galaxy_symmetrised/ED_corr_00_splitbyspatial_central.txt'%snapshot).T
ed_sym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/symmetrised/ED_corr_00.txt'%snapshot).T
ed_msym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/altcent/symmetrised/ED_corr_00_splitbymost_massive.txt'%snapshot).T



ee=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/fiducial/EE_corr_00.txt'%snapshot).T
ee12=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/fiducial/12bins/EE_corr_00.txt'%snapshot).T
ee12_sym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/symmetrised/12bins/EE_corr_00.txt'%snapshot).T
ee_gsym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/galaxy_symmetrised/EE_corr_00_splitbyspatial_central.txt'%snapshot).T
ee_sym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/symmetrised/EE_corr_00.txt'%snapshot).T
ee_red=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/reduced/EE_corr_00_splitbyspatial_central.txt'%snapshot).T
ee_red_sym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/reduced/symmetrised/EE_corr_00_splitbyspatial_central.txt'%snapshot).T
ee_msym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/altcent/symmetrised/EE_corr_00_splitbymost_massive.txt'%snapshot).T

gi=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/fiducial/GIplus_proj_corr_00.txt'%snapshot).T
gi12=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/fiducial/12bins/GIplus_proj_corr_00.txt'%snapshot).T
gi12_sym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/symmetrised/12bins/GIplus_proj_corr_00.txt'%snapshot).T
gi_sym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/symmetrised/GIplus_proj_corr_00.txt'%snapshot).T
gi_gsym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/galaxy_symmetrised/GIplus_proj_corr_00_splitbyspatial_central.txt'%snapshot).T
gi_red=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/reduced/GIplus_proj_corr_00_splitbyspatial_central.txt'%snapshot).T
gi_red_sym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/reduced/symmetrised/GIplus_proj_corr_00_splitbyspatial_central.txt'%snapshot).T
gi_msym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/altcent/symmetrised/GIplus_proj_corr_00_splitbymost_massive.txt'%snapshot).T

ii=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/fiducial/IIplus_proj_corr_00.txt'%snapshot).T
ii12=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/fiducial/12bins/IIplus_proj_corr_00.txt'%snapshot).T
ii12_sym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/symmetrised/12bins/IIplus_proj_corr_00.txt'%snapshot).T
ii_sym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/symmetrised/IIplus_proj_corr_00.txt'%snapshot).T
ii_gsym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/galaxy_symmetrised/IIplus_proj_corr_00_splitbyspatial_central.txt'%snapshot).T
ii_red=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/reduced/IIplus_proj_corr_00_splitbyspatial_central.txt'%snapshot).T
ii_red_sym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/reduced/symmetrised/IIplus_proj_corr_00_splitbyspatial_central.txt'%snapshot).T
ii_msym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot%d/altcent/symmetrised/IIplus_proj_corr_00_splitbymost_massive.txt'%snapshot).T

fig=plt.figure()
ax=fig.add_subplot(411)
delta = 10**-0.06

plt.errorbar(ed[0]*delta, ed[1], yerr=ed[2], color='purple', marker='D', mec='purple', linestyle='none', capsize=0.2, markersize=3.5, label='Basic (unsymmetrised)')
plt.errorbar(ed_sym[0], ed_sym[1], yerr=ed_sym[2], color='royalblue', mec='royalblue', marker='D', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Basic (symmetrised)')
plt.errorbar(ed_red[0], ed_red[1], yerr=ed_red[2], color='plum', marker='^', linestyle='none', mec='plum', capsize=0.2, markersize=3.5, label='Reduced (unsymmetrised)')
plt.errorbar(ed_red_sym[0]/delta, ed_red_sym[1], yerr=ed_red_sym[2], color='hotpink', mec='hotpink', marker='^', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Reduced (symmetrised)')
ax.set_yscale('log', nonposy="clip")
plt.legend(loc='lower left', fontsize=6)
plt.xlim(0.1,33.)
plt.ylim(1e-5,1e-1)
plt.yticks([1e-4, 1e-3,1e-2,1e-1], fontsize=13)
plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=14)
plt.ylabel('ED', fontsize=14)
plt.xticks(visible=False)
plt.xscale('log')

ax=fig.add_subplot(412)

plt.errorbar(ee[0]*delta, ee[1], yerr=ee[2], color='purple', marker='D', mec='purple', linestyle='none', capsize=0.2, markersize=3.5, label='Basic (unsymmetrised)')
plt.errorbar(ee_sym[0], ee_sym[1], yerr=ee_sym[2], color='royalblue', mec='royalblue', marker='D', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Basic (symmetrised)')
plt.errorbar(ee_red[0], ee_red[1], yerr=ee_red[2], color='plum', marker='^', linestyle='none', mec='plum', capsize=0.2, markersize=3.5, label='Reduced (unsymmetrised)')
plt.errorbar(ee_red_sym[0]/delta, ee_red_sym[1], yerr=ee_red_sym[2], color='hotpink', mec='hotpink', marker='^', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Reduced (symmetrised)')
ax.set_yscale('log', nonposy="clip")
plt.xlim(0.1,33.)
plt.ylim(1e-6,3e-1)
plt.yticks([1e-5,1e-3,1e-1], fontsize=13)
plt.xticks(visible=False)
plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=14)
plt.ylabel('EE', fontsize=14)
plt.xscale('log')

ax=fig.add_subplot(413)
plt.errorbar(gi[0]*delta, gi[1], yerr=gi[2], color='purple', marker='D', mec='purple', linestyle='none', capsize=0.2, markersize=3.5, label='Basic (unsymmetrised)')
plt.errorbar(gi_sym[0], gi_sym[1], yerr=gi_sym[2], color='royalblue', mec='royalblue', marker='D', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Basic (symmetrised)')
plt.errorbar(gi_red[0], gi_red[1], yerr=gi_red[2], color='plum', marker='^', linestyle='none', mec='plum', capsize=0.2, markersize=3.5, label='Reduced (unsymmetrised)')
plt.errorbar(gi_red_sym[0]/delta, gi_red_sym[1], yerr=gi_red_sym[2], color='hotpink', mec='hotpink', marker='^', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Reduced (symmetrised)')
ax.set_yscale('log', nonposy="clip")
plt.xlim(0.1,33.)
plt.ylim(5e-6,2e2)
plt.yticks([1e-4, 1e-2,1e0,1e2], fontsize=13)
plt.ylabel(r'$w_{g+}$', fontsize=14)
plt.xticks(visible=False)
plt.xscale('log')

ax=fig.add_subplot(414)
plt.errorbar(ii[0]*delta, ii[1], yerr=ii[2], color='purple', marker='D', mec='purple', linestyle='none', capsize=0.2, markersize=3.5, label='Basic (unsymmetrised)')
plt.errorbar(ii_sym[0], ii_sym[1], yerr=ii_sym[2], color='royalblue', mec='royalblue', marker='D', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Basic (symmetrised)')
plt.errorbar(ii_red[0], ii_red[1], yerr=ii_red[2], color='plum', marker='^', linestyle='none', mec='plum', capsize=0.2, markersize=3.5, label='Reduced (unsymmetrised)')
plt.errorbar(ii_red_sym[0]/delta, ii_red_sym[1], yerr=ii_red_sym[2], color='hotpink', mec='hotpink', marker='^', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Reduced (symmetrised)')
ax.set_yscale('log', nonposy="clip")
plt.xlim(0.1,33.)
plt.ylim(5e-7,9)
plt.yticks([1e-6, 1e-4, 1e-2,1e0], fontsize=13)
plt.ylabel(r'$w_{++}$', fontsize=14)
plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=14)
plt.xscale('log')



plt.subplots_adjust(hspace=0,wspace=0, bottom=0.14, left=0.14, right=0.5, top=0.97)
plt.savefig('/Users/hattifattener/Desktop/corrs-all-mbii-symmetrised-split-reduced.pdf')
xf = np.linspace(ee[0].min(), ee[0].max(), 100)

mask=abs(ee12)[1]>1e-5
p_ee,c_ee=opt.curve_fit(func, ee12[0][mask], ee12[1][mask], sigma=ee12[2][mask])
ee_smooth = func(ee[0], p_ee[0], p_ee[1])
ee_smooth12 = func(ee12[0], p_ee[0], p_ee[1])
p_ed,c_ed=opt.curve_fit(func, ed12[0], ed12[1], sigma=ed12[2])
ed_smooth = func(ed[0], p_ed[0], p_ed[1])
p_gi,c_gi=opt.curve_fit(func, gi12[0], gi12[1], sigma=gi12[2])
gi_smooth = func(gi[0], p_gi[0], p_gi[1])
p_ii,c_ii=opt.curve_fit(func, ii12[0], ii12[1], sigma=ii12[2])
ii_smooth = func(ii[0], p_ii[0], p_ii[1])

p_ed,c_ed=opt.curve_fit(func, ed12[0], ed12[1], sigma=ed12[2])
ed_smooth12 = func(ed12[0], p_ed[0], p_ed[1])



p_gi,c_gi=opt.curve_fit(func, gi12[0], gi12[1], sigma=gi12[2])
gi_smooth12 = func(gi12[0], p_gi[0], p_gi[1])

p_ii,c_ii=opt.curve_fit(func, ii12[0], ii12[1], sigma=ii12[2])
ii_smooth12 = func(ii12[0], p_ii[0], p_ii[1])



mask = (ee[1]-ee_sym[1])>0
p_dee,c_dee = opt.curve_fit(func, ee[0][mask], ee[1][mask]-ee_sym[1][mask], sigma=abs(ee[1][mask]-ee_sym[1][mask])*ee[2][mask]/ee[1][mask])
dee_smooth = func(xf, p_dee[0], p_dee[1])
p_ded,c_ded = opt.curve_fit(func, ed[0], ed[1]-ed_sym[1], sigma=abs(ed[1]-ed_sym[1])*ed[2]/ed[1])
ded_smooth = func(xf, p_ded[0], p_ded[1])
p_dgi,c_dgi = opt.curve_fit(func, gi[0], gi[1]-gi_sym[1], sigma=abs(gi[1]-gi_sym[1])*gi[2]/gi[1])
dgi_smooth = func(xf, p_dgi[0], p_dgi[1])
p_dii,c_dii = opt.curve_fit(func, ii[0], ii[1]-ii_sym[1], sigma=abs(ii[1]-ii_sym[1])*ii[2]/ii[1])
dii_smooth = func(xf, p_dii[0], p_dii[1])


#Reduced
p_ee_red,c_ee_red=opt.curve_fit(func, ee_red[0], ee_red[1], sigma=ee_red[2])
ee_smooth_red = func(xf, p_ee_red[0], p_ee_red[1])
p_ed_red,c_ed_red=opt.curve_fit(func, ed_red[0], ed_red[1], sigma=ed_red[2])
ed_smooth_red = func(xf, p_ed_red[0], p_ed_red[1])
p_gi_red,c_gi_red=opt.curve_fit(func, gi_red[0], gi_red[1], sigma=gi_red[2])
gi_smooth_red = func(xf, p_gi_red[0], p_gi_red[1])
p_ii_red,c_ii_red=opt.curve_fit(func, ii_red[0], ii_red[1], sigma=ii_red[2])
ii_smooth_red = func(xf, p_ii_red[0], p_ii_red[1])

mask = (ee_red[1]-ee_red_sym[1])>0
p_dee_red,c_dee_red = opt.curve_fit(func, ee_red[0][mask], ee_red[1][mask]-ee_red_sym[1][mask], sigma=ee_red[2][mask])
dee_smooth_red = func(xf, p_dee_red[0], p_dee_red[1])
p_ded_red,c_ded_red = opt.curve_fit(func, ed_red[0], abs(ed_red[1]-ed_red_sym[1]), sigma=ed_red[2])
ded_smooth_red = func(xf, p_ded_red[0], p_ded_red[1])
p_dgi_red,c_dgi_red = opt.curve_fit(func, gi_red[0], abs(gi_red[1]-gi_red_sym[1]), sigma=gi_red[2])
dgi_smooth_red = func(xf, p_dgi_red[0], p_dgi_red[1])
p_dii_red,c_dii_red = opt.curve_fit(func, ii_red[0], abs(ii_red[1]-ii_red_sym[1]), sigma=ii_red[2])
dii_smooth_red = func(xf, p_dii_red[0], p_dii_red[1])

# galaxy symmetrised

p_dee_red12,c_dee_red12 = opt.curve_fit(func, ee12[0], ee12[1]-ee_gsym[1], sigma=ee12[2])
dee_smooth_gsym = func(xf, p_dee_red12[0], p_dee_red12[1])
p_ded_gsym,c_ded_gsym = opt.curve_fit(func, ed12[0], ed12[1]-ed_gsym[1], sigma=ed12[2])
ded_smooth_gsym = func(xf, p_ded_gsym[0], p_ded_gsym[1])
#p_dgi,c_dgi = opt.curve_fit(func, gi12[0], gi12[1]-gi_gsym[1], sigma=gi12[2])
#dgi_smooth_gsym = func(xf, p_dgi[0], p_dgi[1])
#p_dii,c_dii = opt.curve_fit(func, ii12[0], ii12[1]-ii_gsym[1], sigma=ii12[2])
#dii_smooth_gsym = func(xf, p_dii[0], p_dii[1])

#
#plt.close()
#plt.errorbar(ee[0], ee[1], yerr=ee[2], color='purple', linestyle='none', marker='x')
#plt.plot(ee12[0], ee12[1]-ee_gsym[1], '.', color='royalblue')
#plt.plot(xf, ee_smooth,color='purple')
#plt.plot(xf, dee_smooth_gsym,color='royalblue')
#plt.xscale('log')
#plt.yscale('log')
#plt.savefig('/Users/hattifattener/Documents/ias/mbii/draft/figures/ee-linefit-reduced.png')




plt.close()
plt.subplot(411)
df_red = (ed_red[1]-ed_red_sym[1]) - (ed[1]-ed_sym[1])
df_msym = (ed[1]-ed_msym[1]) - (ed[1]-ed_sym[1])
#abs(ded_smooth_red/ed_smooth_red) - abs(ded_smooth/ed_smooth)
df_gsym = (ed12[1]-ed_gsym[1]) - (ed12[1]-ed12_sym[1])
#abs(ded_smooth_gsym/ed_smooth) - abs(ded_smooth/ed_smooth)
plt.plot(ed_red[0], abs(df_red), color='purple', lw=2, ls='-', label='Reduced Inertia Tensor')
plt.plot(ed12[0], abs(df_gsym), color='pink', lw=2, ls='--', label='Galaxy-Galaxy Symmetrisation')
plt.plot(ed_msym[0], abs(df_msym), color='plum', lw=2, ls=':', label='Alternative Central Definition')
#plt.plot(ed[0]*delta, , color='k', ls=':', label='ED Reduced')
plt.yscale('log')
plt.xlim(0.1,33.)
#plt.ylim(-0.049,0.01)
#plt.yticks([1e-6, 1e-4, 1e-2,1e0], fontsize=13)
plt.ylabel(r'$\delta \Delta \mathrm{ED}$', fontsize=16)
plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=14)
plt.xscale('log')
#plt.axhline(0,color='k', ls=':')
plt.xticks(visible=False)
plt.yticks([1e-4, 1e-2, 1e0], fontsize=12)
plt.ylim(1e-5,1.2e0)

plt.subplot(412)
df_red = (ee_red[1]-ee_red_sym[1]) - (ee[1]-ee_sym[1])
print df_red
#abs(ded_smooth_red/ed_smooth_red) - abs(ded_smooth/ed_smooth)
df_gsym = (ee12[1]-ee_gsym[1]) - (ee12[1]-ee12_sym[1])
df_msym = (ee[1]-ee_msym[1]) - (ee[1]-ee_sym[1])
#abs(ded_smooth_gsym/ed_smooth) - abs(ded_smooth/ed_smooth)
plt.plot(ee_red[0], abs(df_red), color='purple', lw=2, ls='-', label='Reduced Inertia Tensor')
plt.plot(ee12[0], abs(df_gsym), color='pink', lw=2, ls='--', label='Galaxy-Galaxy Symmetrisation')
plt.plot(ee_msym[0], abs(df_msym), color='plum', lw=2, ls=':', label='Alternative Central Definition')
#plt.plot(ed[0]*delta, , color='k', ls=':', label='ED Reduced')
plt.yscale('log')
plt.xlim(0.1,33.)
#plt.ylim(-0.049,0.01)
#plt.yticks([1e-6, 1e-4, 1e-2,1e0], fontsize=13)
plt.ylabel(r'$\delta \Delta \mathrm{EE}$', fontsize=16)
plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=16)
plt.xscale('log')
#plt.axhline(0,color='k', ls=':')
plt.xticks(visible=False)
plt.yticks([1e-4, 1e-2, 1e0], fontsize=12)
plt.legend(loc='upper right', fontsize=8)
plt.ylim(1e-5,1.2e0)

plt.subplot(413)
df_red = (gi_red[1]-gi_red_sym[1]) - (gi[1]-gi_sym[1])
df_msym = (gi[1]-gi_msym[1]) - (gi[1]-gi_sym[1])
#abs(ded_smooth_red/ed_smooth_red) - abs(ded_smooth/ed_smooth)
df_gsym = (gi12[1]-gi_gsym[1]) - (gi12[1]-gi12_sym[1])
#abs(ded_smooth_gsym/ed_smooth) - abs(ded_smooth/ed_smooth)
plt.plot(gi_red[0], abs(df_red), color='purple', lw=2, ls='-', label='Reduced Inertia Tensor')
plt.plot(gi12[0], abs(df_gsym), color='pink', lw=2, ls='--', label='Galaxy-Galaxy Symmetrisation')
plt.plot(gi_msym[0], abs(df_msym), color='plum', lw=2, ls=':', label='Alternative Central Definition')
#plt.plot(ed[0]*delta, , color='k', ls=':', label='ED Reduced')
plt.yscale('log')
plt.xlim(0.1,33.)
#plt.ylim(-0.049,0.01)
#plt.yticks([1e-6, 1e-4, 1e-2,1e0], fontsize=13)
plt.ylabel(r'$\delta \Delta w_{g+}$', fontsize=16)

plt.xscale('log')
#plt.axhline(0,color='k', ls=':')
plt.xticks(visible=False)
plt.yticks([1e-4, 1e-2, 1e0], fontsize=12)
plt.ylim(1e-5,1.2e0)

plt.subplot(414)
df_red = (ii_red[1]-ii_red_sym[1]) - (ii[1]-ii_sym[1])
#abs(ded_smooth_red/ed_smooth_red) - abs(ded_smooth/ed_smooth)
df_gsym = (ii12[1]-ii_gsym[1]) - (ii12[1]-ii12_sym[1])
df_msym = (ii[1]-ii_msym[1]) - (ii[1]-ii_sym[1])
#abs(ded_smooth_gsym/ed_smooth) - abs(ded_smooth/ed_smooth)
plt.plot(ii_red[0], abs(df_red), color='purple', lw=2, ls='-')
plt.plot(ii12[0], abs(df_gsym), color='pink', lw=2, ls='--')
plt.plot(ii_msym[0], abs(df_msym), color='plum', lw=2, ls=':')
#plt.plot(ed[0]*delta, , color='k', ls=':', label='ED Reduced')
plt.yscale('log')
plt.xlim(0.1,33.)
#plt.ylim(-0.049,0.01)
#plt.yticks([1e-6, 1e-4, 1e-2,1e0], fontsize=13)
plt.ylabel(r'$\delta \Delta w_{++}$', fontsize=16)
plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=16)
plt.xscale('log')
#plt.axhline(0,color='k', ls=':')
plt.xticks(fontsize=12)
plt.yticks([1e-4, 1e-2, 1e0], fontsize=12)
plt.ylim(1e-5,1.2e0)

plt.subplots_adjust(hspace=0,wspace=0, bottom=0.145, left=0.15, right=0.65, top=0.97)
plt.savefig('/Users/hattifattener/Documents/ias/mbii/draft/figures/dcorr-vs-r-robustness.pdf')



plt.close()
plt.subplot(411)
df_red = (ed_red[1]-ed_red_sym[1]) - (ed[1]-ed_sym[1])
df_msym = (ed[1]-ed_msym[1]) - (ed[1]-ed_sym[1])
#abs(ded_smooth_red/ed_smooth_red) - abs(ded_smooth/ed_smooth)
df_gsym = (ed12[1]-ed_gsym[1]) - (ed12[1]-ed12_sym[1])
#abs(ded_smooth_gsym/ed_smooth) - abs(ded_smooth/ed_smooth)
plt.plot(ed_red[0], (df_red)/ed_smooth, color='purple', lw=2, ls='-', label='Reduced Inertia Tensor')
plt.plot(ed12[0], (df_gsym)/ed_smooth12, color='pink', lw=2, ls='--', label='Galaxy-Galaxy Symmetrisation')
plt.plot(ed_msym[0], (df_msym)/ed_smooth, color='plum', lw=2, ls=':', label='Alternative Central Definition')
#plt.plot(ed[0]*delta, , color='k', ls=':', label='ED Reduced')
#plt.yscale('log')
plt.xlim(0.1,33.)
#plt.ylim(-0.049,0.01)
#plt.yticks([1e-6, 1e-4, 1e-2,1e0], fontsize=13)
plt.ylabel(r'$\delta f_\mathrm{ED}$', fontsize=16)
plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=14)
plt.xscale('log')
plt.axhline(0,color='k', ls=':')
plt.xticks(visible=False)
plt.legend(loc='lower right', fontsize=8)
plt.yticks([-0.6,-0.4,-0.2,0.], fontsize=12)
plt.ylim(-0.61,0.18)

plt.subplot(412)
df_red = (ee_red[1]-ee_red_sym[1]) - (ee[1]-ee_sym[1])
df_msym = (ee[1]-ee_msym[1]) - (ee[1]-ee_sym[1])
#abs(ded_smooth_red/ed_smooth_red) - abs(ded_smooth/ed_smooth)
df_gsym = (ee12[1]-ee_gsym[1]) - (ee12[1]-ee12_sym[1])
#abs(ded_smooth_gsym/ed_smooth) - abs(ded_smooth/ed_smooth)
plt.plot(ee_red[0], (df_red)/ee_smooth, color='purple', lw=2, ls='-', label='Reduced Inertia Tensor')
plt.plot(ee12[0], (df_gsym)/ee_smooth12, color='pink', lw=2, ls='--', label='Galaxy-Galaxy Symmetrisation')
plt.plot(ee_msym[0], (df_msym)/ee_smooth, color='plum', lw=2, ls=':', label='Alternative Central Definition')
#plt.plot(ed[0]*delta, , color='k', ls=':', label='ED Reduced')
#plt.yscale('log')
plt.xlim(0.1,33.)
#plt.ylim(-0.049,0.01)
#plt.yticks([1e-6, 1e-4, 1e-2,1e0], fontsize=13)
plt.ylabel(r'$\delta f_\mathrm{EE}$', fontsize=16)
plt.xscale('log')
plt.axhline(0,color='k', ls=':')
plt.xticks(visible=False)
plt.yticks([-0.6,-0.4,-0.2,0.], fontsize=12)
plt.ylim(-0.61,0.18)

plt.subplot(413)
df_red = (gi_red[1]-gi_red_sym[1]) - (gi[1]-gi_sym[1])
df_msym = (gi[1]-gi_msym[1]) - (gi[1]-gi_sym[1])
#abs(ded_smooth_red/ed_smooth_red) - abs(ded_smooth/ed_smooth)
df_gsym = (gi12[1]-gi_gsym[1]) - (gi12[1]-gi12_sym[1])
#abs(ded_smooth_gsym/ed_smooth) - abs(ded_smooth/ed_smooth)
plt.plot(gi_red[0], (df_red)/gi_smooth, color='purple', lw=2, ls='-', label='Reduced Inertia Tensor')
plt.plot(gi12[0], (df_gsym)/gi_smooth12, color='pink', lw=2, ls='--', label='Galaxy-Galaxy Symmetrisation')
plt.plot(gi_msym[0], (df_msym)/gi_smooth, color='plum', lw=2, ls=':', label='Alternative Central Definition')
#plt.plot(ed[0]*delta, , color='k', ls=':', label='ED Reduced')
#plt.yscale('log')
plt.xlim(0.1,33.)
#plt.ylim(-0.049,0.01)
#plt.yticks([1e-6, 1e-4, 1e-2,1e0], fontsize=13)
plt.ylabel(r'$\delta f_{g+}$', fontsize=16)
plt.xscale('log')
plt.axhline(0,color='k', ls=':')
plt.xticks(visible=False)
plt.yticks([-0.6,-0.4,-0.2,0.], fontsize=12)
plt.ylim(-0.61,0.18)

plt.subplot(414)
df_red = (ii_red[1]-ii_red_sym[1]) - (ii[1]-ii_sym[1])
df_msym = (ii[1]-ii_msym[1]) - (ii[1]-ii_sym[1])
#abs(ded_smooth_red/ed_smooth_red) - abs(ded_smooth/ed_smooth)
df_gsym = (ii12[1]-ii_gsym[1]) - (ii12[1]-ii12_sym[1])
#abs(ded_smooth_gsym/ed_smooth) - abs(ded_smooth/ed_smooth)
plt.plot(ii_red[0], (df_red)/ii_smooth, color='purple', lw=2, ls='-', label='Reduced Inertia Tensor')
plt.plot(ii12[0], (df_gsym)/ii_smooth12, color='pink', lw=2, ls='--', label='Galaxy-Galaxy Symmetrisation')
plt.plot(ii_msym[0], (df_msym)/ii_smooth, color='plum', lw=2, ls=':', label='Alternative Central Definition')
#plt.plot(ed[0]*delta, , color='k', ls=':', label='ED Reduced')
#plt.yscale('log')
plt.xlim(0.1,33.)
#plt.ylim(-0.049,0.01)
#plt.yticks([1e-6, 1e-4, 1e-2,1e0], fontsize=13)
plt.ylabel(r'$\delta f_\mathrm{g+}$', fontsize=16)
plt.xscale('log')
plt.axhline(0,color='k', ls=':')
plt.xticks(visible=False)
plt.ylim(-0.61,0.18)
plt.ylabel(r'$\delta f_{++}$', fontsize=16)
plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=16)
#plt.axhline(0,color='k', ls=':')
plt.xticks(fontsize=12)
plt.yticks([-0.6,-0.4,-0.2,0.], fontsize=12)
#plt.ylim(1e-5,1.2e0)

plt.subplots_adjust(hspace=0,wspace=0, bottom=0.145, left=0.15, right=0.65, top=0.97)
plt.savefig('/Users/hattifattener/Documents/ias/mbii/draft/figures/df-vs-r-robustness.pdf')
