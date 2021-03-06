import fitsio as fi
import treecorr
import numpy as np 
import argparse
import yaml
#import mbii.lego_tools as util
from mbii.pipeline.twopoint.jackknife import multipurpose_errors as errors 

period={'massiveblackii':100, 'illustris':75}

def compute(options, binning):
	print('Shape data : %s'%options['2pt']['shapes'])

	data = fi.FITS(options['2pt']['shapes'])[-1].read()
	data_sym = fi.FITS(options['2pt']['shapes_symmetrised'])[-1].read()

	if 'npart_cut' in options['2pt'].keys():
		nlow = options['2pt']['npart_cut']
		print('Imposing additional cut at npart>%d'%nlow)
		data = data[(data['npart_dm']>nlow)]
		data_sym = data_sym[(data_sym['npart_dm']>nlow)]
	else:
		nlow=-1

	if 'npart_cut_upper' in options['2pt'].keys():
		nhigh = options['2pt']['npart_cut_upper']
		print('Imposing additional cut at npart<%d'%nhigh)
		data = data[(data['npart_dm']<nhigh)]
		data_sym = data_sym[(data_sym['npart_dm']<nhigh)]
	else:
		nhigh = -1

	splitflag = options['2pt']['split']

	if splitflag is not None:
		name = options['2pt']['split']
		print('Dividing catalogue by %s'%name)
		mask = (data[name]>=options['2pt']['split_val'])
		print('%3.3f percent split'%(data[name][mask].size*100./data[name].size))
	else:
		print('No catalogue split required.')
		mask = np.ones(data.size).astype(bool)

	if options['2pt']['binning']=='log':
		rbins = np.logspace(np.log10(options['2pt']['rmin']), np.log10(options['2pt']['rmax']), binning+1)
	elif options['2pt']['binning']=='equal':
		rbins = util.equalise_binning(data[mask], data[mask], options['2pt']['rmin'], options['2pt']['rmax'], binning+1)

	print('Setting up correlations')

	if splitflag:
		suffix = '_splitby%s'%options['2pt']['split']
	else:
		suffix=''

	if (nlow>=0):
		suffix+='-ndm_part_low%d'%nlow
	if (nhigh>0):
		suffix+='-ndm_part_high%d'%nhigh

	# don't need randoms here if we know the period of the box
	print('Computing correlation functions.')

	if splitflag:
		cat1 = data[mask]
		cat1_sym = data_sym[mask]
		cat2 = data[np.invert(mask)]
		cat2_sym = data_sym[np.invert(mask)]

		print('11')
		F11,R11 = errors.jackknife('ed', cat1, cat1, cat1_sym, cat1_sym, options, nbins=binning)
		export_array('%s/ED_var_11%s.txt'%(options['2pt']['savedir'], suffix), rbins, F11, R11)
		
		print('22')
		F22,R22 = errors.jackknife('ed', cat2, cat2, cat2_sym, cat2_sym, options, nbins=binning)
		export_array('%s/ED_var_22%s.txt'%(options['2pt']['savedir'], suffix), rbins, F22, R22)
	
		print('12')
		F12,R12 = errors.jackknife('ed', cat1, cat2, cat1_sym, cat2_sym, options, nbins=binning)
		print('21')
		F21,R21 = errors.jackknife('ed', cat2, cat1, cat2_sym, cat1_sym, options, nbins=binning)

		export_array('%s/ED_var_12%s.txt'%(options['2pt']['savedir'], suffix), rbins, F12, R12)
		export_array('%s/ED_var_21%s.txt'%(options['2pt']['savedir'], suffix), rbins, F21, R21)

	print('00')
	cat0 = data
	cat0_sym = data_sym
	F00,R00 = errors.jackknife('ed', cat0, cat0, cat0_sym, cat0_sym, options, nbins=binning)
	export_array('%s/ED_var_00%s.txt'%(options['2pt']['savedir'], suffix), rbins, F00, R00)

	print('Done')


def export_array(path, edges, result, error):
	x = np.sqrt(edges[:-1]*edges[1:])
	out = np.vstack((x, result, error))
	print('Exporting to %s'%path)
	np.savetxt(path, out.T)

