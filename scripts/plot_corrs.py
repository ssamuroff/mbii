import pylab as plt
plt.switch_backend('agg')

plt.close()
plt.plot(res[0], res[1][0], color='purple', ls='--', lw=2, label='$ss$ 1h')
plt.plot(res[0], res[1][1], color='purple', ls=':', lw=2, label='$ss$ 2h')
plt.plot(res[0], res[1][2], color='royalblue', ls='--', lw=2, label='$sc$ 1h')
plt.plot(res[0], res[1][3], color='royalblue', ls=':', lw=2, label='$sc$ 2h')
#plt.plot(res[0], res[1][4], color='forestgreen', ls='--', lw=2, label='$cc$ 1h')
plt.plot(res[0], res[1][5], color='forestgreen', ls=':', lw=2, label='$cc$ 2h')

plt.plot(res_nosep[0], res_nosep[1], color='k', ls='-', lw=2)



plt.xscale('log')
plt.yscale('symlog')
plt.axhline(0,color='k', ls='-')
leg = plt.legend(loc='center left')
leg.get_frame().set_alpha(0.6)
plt.xlim(30,5e4)
plt.xlabel('Physical Scale $r$ / $h^{-1}$ kpc', fontsize=18)
plt.ylabel(r'$\xi_{gg}(r)$', fontsize=18)
plt.subplots_adjust(bottom=0.12, top=0.98)
plt.savefig('/home/ssamurof/tmp3.png')
