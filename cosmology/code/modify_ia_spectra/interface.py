from cosmosis.datablock import option_section, names
import scipy.interpolate as spi
import numpy as np
import pdb


def setup(options):
	gipath = options.get_string(option_section, 'gifile', default='')
	iipath = options.get_string(option_section, 'iifile', default='')
	gipath_ref = options.get_string(option_section, 'gireffile', default='')
	iipath_ref = options.get_string(option_section, 'iireffile', default='')

	gi,ii,gi_ref,ii_ref='','','',''

	if len(gipath)>0:
		gi = np.loadtxt(gipath)
	if len(iipath)>0:
		gi_ref = np.loadtxt(gipath_ref)
	if len(gipath_ref)>0:
		ii = np.loadtxt(iipath)
	if len(iipath_ref)>0:
		ii_ref = np.loadtxt(iipath_ref)

	config = gi, gi_ref, ii, ii_ref
	return config

def execute(block, config):
	Pgi, Pgi0, Pii, Pii0 = config

	gi_interpolator = spi.interp1d(Pgi[0], Pgi[1]) 
	gi_interpolator0 = spi.interp1d(Pgi0[0], Pgi0[1]) 

	ii_interpolator = spi.interp1d(Pii[0], Pii[1]) 
	ii_interpolator0 = spi.interp1d(Pii0[0], Pii0[1])

	GI = block['matter_intrinsic_power', 'p_k']
	II = block['intrinsic_power', 'p_k']
	k = block['matter_intrinsic_power', 'k_h']

	dgi = (gi_interpolator(k)-gi_interpolator0(k))/gi_interpolator0(k)
	dii = (ii_interpolator(k)-ii_interpolator0(k))/ii_interpolator0(k)

	mask = (k<1) & (k>=0.2)
	dgi[mask]*=np.exp(k[mask]*4)/np.exp(k[mask]*4)[-1]
	dii[mask]*=np.exp(k[mask]*4)/np.exp(k[mask]*4)[-1]
	dgi[(k<0.2)]*=0
	dii[(k<0.2)]*=0

	GI_modified = np.array([spect*(1-dgi) for spect in GI])
	II_modified = np.array([spect*(1-dii) for spect in II])

	block['matter_intrinsic_power', 'p_k'] = GI_modified
	block['intrinsic_power', 'p_k'] = II_modified

	import pdb ; pdb.set_trace()

	print 'Overwriting IA power spectra.'

	return 0

def cleanup(config):
	return 0
