import pylab as plt
import numpy as np
import fitsio as fi
import argparse

plt.switch_backend('pdf')


parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--verbosity', '-v', type=int, action='store', default=1)
parser.add_argument('--catalogue', '-c', type=str, action='store')
args = parser.parse_args()

data = fi.FITS(args.catalogue)[-1].read()


# Basic mass distributions, split into centrals and satellites
Mm = data['matter_mass']
Mb = data['baryon_mass']
csmask = data['central']
plt.close()
H,b=np.histogram(np.log10(Mm), bins=np.linspace(-0.3,2.5,45), normed=1)
Hc,bc=np.histogram(np.log10(Mm[(csmask==1)]), bins=np.linspace(-0.3,2.5,45), normed=1)
Hs,bs=np.histogram(np.log10(Mm[(csmask!=1)]), bins=np.linspace(-0.3,2.5,45), normed=1)
xb=(b[1:]+b[:-1])/2
plt.plot(xb,H, color='purple', lw=2, label='All Galaxies')
plt.fill_between(xb,Hs, color='royalblue', alpha=0.2, lw=2, label='Satellites')
plt.fill_between(xb,Hc, color='red', alpha=0.2, lw=2, label='Centrals')
plt.legend(loc='upper right')
plt.ylim(ymin=0)
plt.xlabel('Subhalo Matter Mass $\log(M_\mathrm{m})$ / $10^{10}h^{-1} M_{*}$')
plt.subplots_adjust(bottom=0.14, top=0.98, right=0.98, left=0.14)
plt.savefig('/physics2/ssamurof/massive_black_ii/plots/diagnostics/subhalo_mmass-cs_decomp-v3.png')

plt.close()
H,b=np.histogram(np.log10(Mb), bins=np.linspace(-1.7,2.5,50), normed=1)
Hc,bc=np.histogram(np.log10(Mb[(csmask==1)]), bins=np.linspace(-1.7,2.5,50), normed=1)
Hs,bs=np.histogram(np.log10(Mb[(csmask!=1)]), bins=np.linspace(-1.7,2.5,50), normed=1)
xb=(b[1:]+b[:-1])/2
plt.plot(xb,H, color='purple', lw=2, label='All Galaxies')
plt.fill_between(xb,Hs, color='royalblue', alpha=0.2, lw=2, label='Satellites')
plt.fill_between(xb,Hc, color='red', alpha=0.2, lw=2, label='Centrals')
plt.legend(loc='upper right')
plt.ylim(ymin=0)
plt.xlabel('Subhalo Baryonic Mass $\log(M_\mathrm{b})$ / $10^{10} h^{-1} M_{*}$')
plt.subplots_adjust(bottom=0.14, top=0.98, right=0.98, left=0.14)
plt.savefig('/physics2/ssamurof/massive_black_ii/plots/diagnostics/subhalo_bmass-cs_decomp-v3.png')





# Basic mass distributions, split into discs and bulges
Mm = data['matter_mass']
Mb = data['baryon_mass']
R=(data['vrot']/data['sigma'])
bdmask = R<np.median(R)
plt.close()
H,b=np.histogram(np.log10(Mm), bins=np.linspace(-0.3,2.5,45), normed=1)
Hc,bc=np.histogram(np.log10(Mm[(bdmask)]), bins=np.linspace(-0.3,2.5,45), normed=1)
Hs,bs=np.histogram(np.log10(Mm[(bdmask)]), bins=np.linspace(-0.3,2.5,45), normed=1)
xb=(b[1:]+b[:-1])/2
plt.plot(xb,H, color='purple', lw=2, label='All Galaxies')
plt.fill_between(xb,Hs, color='royalblue', alpha=0.2, lw=2, label='Rotation Dominated')
plt.fill_between(xb,Hc, color='red', alpha=0.2, lw=2, label='Pressure Dominated')
plt.legend(loc='upper right')
plt.ylim(ymin=0)
plt.xlabel('Subhalo Matter Mass $\log(M_\mathrm{m})$ / $10^{10}h^{-1} M_{*}$')
plt.subplots_adjust(bottom=0.14, top=0.98, right=0.98, left=0.14)
plt.savefig('/physics2/ssamurof/massive_black_ii/plots/diagnostics/subhalo_mmass-bd_decomp-v3.png')

plt.close()
H,b=np.histogram(np.log10(Mb), bins=np.linspace(-1.7,2.5,50), normed=1)
Hc,bc=np.histogram(np.log10(Mb[(bdmask==1)]), bins=np.linspace(-1.7,2.5,50), normed=1)
Hs,bs=np.histogram(np.log10(Mb[(bdmask!=1)]), bins=np.linspace(-1.7,2.5,50), normed=1)
xb=(b[1:]+b[:-1])/2
plt.plot(xb,H, color='purple', lw=2, label='All Galaxies')
plt.fill_between(xb,Hs, color='royalblue', alpha=0.2, lw=2, label='Rotation Dominated')
plt.fill_between(xb,Hc, color='red', alpha=0.2, lw=2, label='Pressure Dominated')
plt.legend(loc='upper right')
plt.ylim(ymin=0)
plt.xlabel('Subhalo Baryonic Mass $\log(M_\mathrm{b})$ / $10^{10} h^{-1} M_{*}$')
plt.subplots_adjust(bottom=0.14, top=0.98, right=0.98, left=0.14)
plt.savefig('/physics2/ssamurof/massive_black_ii/plots/diagnostics/subhalo_bmass-bd_decomp-v3.png')

