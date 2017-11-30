import mbii.lego_tools as util
import pylab as plt
plt.switch_backend('agg')

do_errorbars=False

plt.close()
plt.plot(R, res[0], color='purple', ls='--', lw=2, label='$ss$ 1h')
B=util.process_bootstrap(bootstrap_iterations[0])
if do_errorbars: plt.fill_between(res[0], res[1][0]-B, res[1][0]+B, color='purple', alpha=0.3)

B=util.process_bootstrap(bootstrap_iterations[1])
plt.plot(R, res[1], color='purple', ls=':', lw=2, label='$ss$ 2h')
if do_errorbars: plt.fill_between(res[0], res[1][1]-B, res[1][1]+B, color='purple', alpha=0.3)

B=util.process_bootstrap(bootstrap_iterations[2])
plt.plot(R, res[2], color='royalblue', ls='--', lw=2, label='$sc$ 1h')
if do_errorbars: plt.fill_between(res[0], res[1][2]-B, res[1][2]+B, color='royalblue', alpha=0.3)

B=util.process_bootstrap(bootstrap_iterations[3])
plt.plot(R, res[3], color='royalblue', ls=':', lw=2, label='$sc$ 2h')
if do_errorbars: plt.fill_between(res[0], res[1][3]-B, res[1][3]+B, color='royalblue', alpha=0.3)

B=util.process_bootstrap(bootstrap_iterations[5])
plt.plot(R, res[5], color='forestgreen', ls=':', lw=2, label='$cc$ 2h')
if do_errorbars: plt.fill_between(res[0], res[1][5]-B, res[1][5]+B, color='forestgreen', alpha=0.3)

plt.plot(res_nosep[0], res_nosep[1], color='k', ls='-', lw=1, alpha=0.5)



plt.xscale('log')
plt.yscale('symlog')
plt.axhline(0,color='k', ls='-')
leg = plt.legend(loc='upper right')
leg.get_frame().set_alpha(0.6)
plt.xlim(30,5e4)
plt.xlabel('Physical Scale $r$ / $h^{-1}$ kpc', fontsize=18)
plt.ylabel(r'$\xi_{gg}(r)$', fontsize=18)
plt.subplots_adjust(bottom=0.12, top=0.98)
plt.savefig('/home/ssamurof/tmp3.png')




import mbii.lego_tools as util
import pylab as plt
plt.switch_backend('agg')

# This just compares the mixed 1h, 2h, satellite, central position correlation function
# as calculated using treecorr and halotools
plt.close('all')
fig = plt.figure()
plt.plot(R_nosep, res_nosep, color='hotpink', ls='-', lw=1, alpha=1, label='halotools')
plt.plot(R_nosep_treecorr, res_nosep_treecorr, color='purple', ls='--', lw=1, alpha=1, label='treecorr')



plt.xscale('log')
plt.yscale('symlog')
plt.axhline(0,color='k', ls='-')
leg = plt.legend(loc='upper right')
leg.get_frame().set_alpha(0.6)
plt.xlim(30,5e4)
#plt.ylim(1e9,5e13)
plt.xlabel('Physical Scale $r$ / $h^{-1}$ kpc', fontsize=18)
plt.ylabel(r'$\xi_{gg}(r)$', fontsize=18)
plt.subplots_adjust(bottom=0.12, top=0.98)
plt.savefig('/home/ssamurof/tmp4.png')
