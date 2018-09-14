import numpy as np
import scipy.optimize as opt
import pylab as plt
plt.switch_backend('agg')
plt.style.use('y1a1')

snapshot = 85

def func(x, a, b):
	return 10**b * x**a

f={}
for corr in ['ee', 'ed', 'gi', 'ii']: f[corr]=[]
nbins = np.array([300,645,1000, 1239,1954,2700,3775,4842,6673,10379])


# Load the data for this snapshot
for nlow in nbins:
	print nlow
	plt.close()
	ed=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/fiducial/ED_corr_00-ndm_part_low%d.txt'%nlow).T
	ed_sym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/symmetrised/ED_corr_00-ndm_part_low%d.txt'%nlow).T
	ee=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/fiducial/EE_corr_00-ndm_part_low%d.txt'%nlow).T
	ee_sym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/symmetrised/EE_corr_00-ndm_part_low%d.txt'%nlow).T
	gi=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/fiducial/GIplus_proj_corr_00-ndm_part_low%d.txt'%nlow).T
	gi_sym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/symmetrised/GIplus_proj_corr_00-ndm_part_low%d.txt'%nlow).T
	ii=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/fiducial/IIplus_proj_corr_00-ndm_part_low%d.txt'%nlow).T
	ii_sym=np.loadtxt('/Users/hattifattener/Documents/ias/mbii/2pt/ns300_nd1000/v10/snapshot85/symmetrised/IIplus_proj_corr_00-ndm_part_low%d.txt'%nlow).T

	fig=plt.figure()
	ax=fig.add_subplot(411)
	delta = 10**-0.06

	plt.errorbar(ed[0]*delta, ed[1], yerr=ed[2], color='purple', marker='D', mec='purple', linestyle='none', capsize=0.2, markersize=3.5, label='Unsymmetrised')
	plt.errorbar(ed_sym[0], ed_sym[1], yerr=ed_sym[2], color='royalblue', mec='royalblue', marker='D', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Symmetrised')
	ax.set_yscale('log', nonposy="clip")
	plt.legend(loc='lower left', fontsize=8)
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
	ax.set_yscale('log', nonposy="clip")
	plt.xlim(0.1,33.)
	plt.ylim(1e-6,3e-1)
	plt.yticks([1e-5,1e-3,1e-1], fontsize=13)
	plt.xticks(visible=False)
	#plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=14)
	plt.ylabel('EE', fontsize=14)
	plt.xscale('log')


	ax=fig.add_subplot(413)
	plt.errorbar(gi[0]*delta, gi[1], yerr=gi[2], color='purple', marker='D', mec='purple', linestyle='none', capsize=0.2, markersize=3.5, label='Basic (unsymmetrised)')
	plt.errorbar(gi_sym[0], gi_sym[1], yerr=gi_sym[2], color='royalblue', mec='royalblue', marker='D', mfc='none', linestyle='none', capsize=0.2, markersize=3.5, label='Basic (symmetrised)')
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
	ax.set_yscale('log', nonposy="clip")
	plt.xlim(0.1,33.)
	plt.ylim(5e-7,9)
	plt.yticks([1e-6, 1e-4, 1e-2,1e0], fontsize=13)
	plt.ylabel(r'$w_{++}$', fontsize=14)
	plt.xlabel('$r$ / $h^{-1}$ Mpc', fontsize=14)
	plt.xscale('log')

	plt.subplots_adjust(hspace=0,wspace=0, right=0.5, top=0.98, bottom=0.14)
	plt.savefig('/Users/hattifattener/Documents/ias/mbii/draft/figures/corrs-all-mbii-symmetrised-split-snapshot%d-nlow%d.pdf'%(snapshot, nlow))

	xf = np.linspace(gi[0].min(), gi[0].max(), 100)

	p_ee,c_ee=opt.curve_fit(func, ee[0], ee[1], sigma=ee[2])
	ee_smooth = func(ee[0], p_ee[0], p_ee[1])
	p_ed,c_ed=opt.curve_fit(func, ed[0], ed[1], sigma=ed[2])
	ed_smooth = func(ed[0], p_ed[0], p_ed[1])
	p_gi,c_gi=opt.curve_fit(func, gi[0], gi[1], sigma=gi[2])
	gi_smooth = func(ed[0], p_gi[0], p_gi[1])
	p_ii,c_ii=opt.curve_fit(func, ii[0], ii[1], sigma=ii[2])
	ii_smooth = func(ed[0], p_ii[0], p_ii[1])
#
	p_dee,c_dee = opt.curve_fit(func, ee[0], ee[1]-ee_sym[1], sigma=ee[2])
	dee_smooth = func(ee[0], p_dee[0], p_dee[1])
	p_ded,c_ded = opt.curve_fit(func, ed[0], ed[1]-ed_sym[1], sigma=ed[2])
	ded_smooth = func(ed[0], p_ded[0], p_ded[1])
	p_dgi,c_dgi = opt.curve_fit(func, gi[0], gi[1]-gi_sym[1], sigma=gi[2])
	dgi_smooth = func(ed[0], p_dgi[0], p_dgi[1])
	p_dii,c_dii = opt.curve_fit(func, ii[0], ii[1]-ii_sym[1], sigma=ii[2])
	dii_smooth = func(ed[0], p_dii[0], p_dii[1])
#

	#if snapshot==73:
	#	import pdb ; pdb.set_trace()

	plt.close()
	plt.subplot(211)
	#p1, = plt.plot(gg0[0], abs(gg0[1]-gg1[1]), '.', color='steelblue', ls=':', label=r'$\xi_{gg}$')
	p2, = plt.plot(ed[0], abs(ed[1]-ed_sym[1]), color='purple', ls='--', label='ED')
	p3, = plt.plot(ee[0], abs(ee[1]-ee_sym[1]), color='royalblue', ls=':', label='EE')
	p4, = plt.plot(gi[0], abs(gi[1]-gi_sym[1]), color='plum', ls='-.', label='$w_{g+}$')
	p5, = plt.plot(ii[0], abs(ii[1]-ii_sym[1]), color='pink', ls='-', label='$w_{++}$')
	plt.xscale('log')
	plt.yscale('log')
	plt.ylim(4e-5, 3e-1)
	plt.xlim(0.1,33)
	plt.ylabel('Absolute Residual', fontsize=16)
	plt.yticks(fontsize=16)
	l1 = plt.legend([p2, p3], ['ED', 'EE'], loc='upper right', fontsize=14)
	l2 = plt.legend([p4, p5], ['$w_{g+}$', '$w_{++}$'], loc='lower left', fontsize=14)
	plt.gca().add_artist(l1)
	#plt.legend(loc='upper right', fontsize=12)
	#plt.axhline(0, color='k', ls=':')
	plt.xticks(visible=False)

	plt.subplot(212)
	#plt.plot(gg0[0], (gg0[1]-gg1[1])/gg1[1], '.', color='steelblue', ls=':', label=r'$\xi_{gg}$')
	
	plt.plot(ed[0], (ed[1]-ed_sym[1])/ed[1], 'D', color='purple', ls='none', label='ED')
	plt.plot(ee[0], (ee[1]-ee_sym[1])/ee[1], 'D', color='royalblue', ls='none', label='EE')
	plt.plot(gi[0], (gi[1]-gi_sym[1])/gi[1], 'D', color='plum', ls='none', label='$w_{g+}$')
	plt.plot(ii[0], (ii[1]-ii_sym[1])/ii[1], 'D', color='pink', ls='none', label='$w_{++}$')
	plt.plot(ed[0], ded_smooth/ed_smooth, color='purple', ls='-.', label='ED')
	plt.plot(ee[0], dee_smooth/ee_smooth, color='royalblue', ls='-.', label='EE')
	plt.plot(gi[0], dgi_smooth/gi_smooth, color='plum', ls='-.', label='$w_{g+}$')
	plt.plot(ii[0], dii_smooth/ii_smooth, color='pink', ls='-', label='$w_{++}$')
	plt.xscale('log')

	for corr in ['ee', 'ed', 'gi', 'ii']: #['ee', 'ed', 'gi', 'ii']: 
		exec('R = ((%s[1]-%s_sym[1])/%s_smooth)'%(corr,corr,corr))
		f[corr].append(R[(gi[0]>2.0)].mean())

	plt.ylim(-0.15, 0.79)
	plt.xlim(0.1,33)
	plt.yticks(fontsize=16)
	plt.ylabel('Fractional Residual', fontsize=16)
	plt.xticks(fontsize=16)
	plt.xlabel('$r$ / $h{-1}$ Mpc', fontsize=18)
	plt.axhline(0, color='k', ls=':')

	plt.subplots_adjust(hspace=0, left=0.14, right=0.7, bottom=0.14, top=0.97)
	plt.savefig('/Users/hattifattener/Documents/ias/mbii/draft/figures/all-corrs-symdiff-snapshot%d-nlow%d.pdf'%(snapshot, nlow))




Fee = f['ee']
Fed = f['ed']
Fgi = f['gi']
Fii = f['ii']

plt.close()

plt.plot(nbins, Fed, marker='.', ms=5.5, color='purple', linestyle='none', label='ED')
plt.plot(nbins+75, Fee, marker='D', ms=5.5, color='steelblue', linestyle='none', label='EE')
plt.plot(nbins-75, Fgi, marker='x', ms=5.5, color='plum', linestyle='none', label='$w_{g+}$')
plt.plot(nbins-75, Fii, marker='^', ms=5.5, color='pink', linestyle='none', label='$w_{++}$')
plt.xlim(150,5000)
plt.legend(loc='upper left')
plt.xticks(fontsize=16)
plt.xlabel(r'Lower Cut $>N_\mathrm{m}$', fontsize=18)
plt.ylabel('Mean Residual $f(z | r>2 h^{-1}$ Mpc$)$', fontsize=16)
plt.axhline(0, color='k', ls=':')
plt.ylim(-0.06, 1.15)
plt.subplots_adjust(hspace=0, left=0.14, right=0.98, bottom=0.14, top=0.97)
plt.savefig('/Users/hattifattener/Documents/ias/mbii/draft/figures/fres-vs-nlow-allcorrs.pdf')

plt.close()

import fitsio as fi
mdm = 1.1e7
mbins = nbins*mdm

cat = fi.FITS('/Users/hattifattener/Documents/ias/mbii/cats/final2/postprocessed/base_subhalo_shapes-v10-ndm300-nst300-snapshot85.fits')[-1].read()
mdata = cat['npart_dm']*mdm

fig, ax1 = plt.subplots()

ax1.plot(mbins/1e10, Fed, marker='.', ms=5.5, color='purple', linestyle='none', label='ED')
plt.plot(mbins/1e10, Fee, marker='D', ms=5.5, color='steelblue', linestyle='none', label='EE')
ax1.plot(mbins/1e10, Fgi, marker='x', ms=5.5, color='plum', linestyle='none', label='$w_{g+}$')
ax1.plot(mbins/1e10, Fii, marker='^', ms=5.5, color='pink', linestyle='none', label='$w_{++}$')
ax1.set_xlim(100*mdm/1e10, 13.99 )
plt.legend(loc='upper right')
plt.xticks(fontsize=16)
plt.xlabel(r'Lower Mass Cut $>M_\mathrm{m} / 10^{10} h^{-1} M_\odot$', fontsize=18)
ax1.set_ylabel('Mean Residual $f(z | r>2 h^{-1}$ Mpc$)$', fontsize=16)
ax1.axhline(0, color='k', ls=':')
ax1.set_ylim(-0.03, 1.17)
plt.axvline(1000*mdm/1e10, color='k', ls=':')

ax2 = ax1.twinx()
m0 = 300 * mdm / 1e10
m1 = 13000 * mdm / 1e10
H,b = np.histogram(mdata/1e10, bins=np.linspace(m0,m1, 180))
x = (b[:-1]+b[1:])/2
ax2.fill_between(x, H, color='purple', alpha=0.2)
ax2.set_ylim(ymin=0, ymax=H.max()+10)
ax2.set_xlim(100*mdm/1e10, 13.99 )
plt.subplots_adjust(hspace=0, left=0.14, right=0.98, bottom=0.148, top=0.97)
plt.savefig('/Users/hattifattener/Documents/ias/mbii/draft/figures/fres-vs-mlow-allcorrs.pdf')