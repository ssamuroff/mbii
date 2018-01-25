import numpy as np
import argparse
import yaml
import fitsio as fi
import scipy.spatial as sps
from mbii.readsubhalo import *
from mbii.properties import *
import pymysql as mdb
import numpy as np
from numpy.core.records import fromarrays


dt = [
      ('object_id', int), # unique ID for each galaxy 
      ('baryon_mass', float), # What it sounds like, in units of 10^10 h^-1 M_*
      ('matter_mass', float),  # What it sounds like, in units of 10^10 h^-1 M_*
      ('npart_dm', int), # Number of DM particles in the subhalo
      ('npart_baryon', int), # Number of star particles in the subhalo
      ('e1', float), ('e2', float), # Projected ellipticities
      ('a1', float), ('a2', float), ('a3', float), # Normalised 3 vector along the major axis 
      ('b1', float), ('b2', float), ('b3', float),
      ('c1', float), ('c2', float), ('c3', float), # Normalised 3 vector along the major axis 
      ('e1_dm', float), ('e2_dm', float), # Projected ellipticities, defined by the dark matter distribution
      ('a1_dm', float), ('a2_dm', float), ('a3_dm', float), # Normalised 3 vector along the major axis, defined by the dark matter distribution 
      ('b1_dm', float), ('b2_dm', float), ('b3_dm', float),
      ('c1_dm', float), ('c2_dm', float), ('c3_dm', float), # Normalised 3 vector along the major axis , defined by the dark matter distribution
      ('x', float), ('y', float), ('z', float), # 3D position in Cartesian coordinates, in Mpc h^-1
      ('central', int), # Flag to indicate central galaxies
      ('rh', int), # Radial distance from the centre of the halo
      ('halo_id', int), # ID of the halo in which each galaxy resides
      ('nocc', int)] # Number of galaxies in that halo (so for any galaxy there are nocc-1 others sharing the same halo)


class halo_wrapper:
	def __init__(self, snapshot, verbosity=1):
		# This is (hopefully) the one and only time we need to call on the coma DB in the catalogue pipeline

		sqlserver='localhost'
		user='flanusse'
		password='mysqlpass'
		dbname='mb2_hydro'
		unix_socket='/home/rmandelb.proj/flanusse/mysql/mysql.sock'
		db = mdb.connect(sqlserver, user, password, dbname, unix_socket=unix_socket)

		# Setup the database connection and query for the centroids of all the groups
		c = db.cursor()
		sql = 'SELECT x,y,z,groupId,mass,len FROM subfind_halos WHERE snapnum=85;'
		if verbosity>0:
			print 'Submitting query...'
		c.execute(sql)
		self.info = fromarrays(np.array(c.fetchall()).squeeze().T,names='x,y,z,groupId,mass,len')

		# Convert to Mpc h^-1
		for name in ['x', 'y', 'z']:
			self.info[name]/=1e3

		self.ids = np.arange(0, len(self.info), 1)

		if verbosity>0:
			print 'Done.'

		self.verbosity = verbosity

	def populate(self, data):
		# Match up each galaxy (subhalo) to the nearest halo centroid
		if self.verbosity>0:
			print 'Building halo KD tree'
		# Build the tree
		xyz0 = np.array([self.info['x'], self.info['y'], self.info['z']])
		tree = sps.KDTree(xyz0.T)

		# Query it for each 3D galaxy position
		if self.verbosity>0:
			print 'Querying tree'
		xyz = np.array([data['x'], data['y'], data['z']])
		R,ind = tree.query(xyz.T, k=1)

		return R, ind, self.info[ind]


class catalogue:
	def __init__(self):
		print 'Set up subhalo catalogue'

	def build_selection_mask(self, cuts, baryons, dm, verbosity=1):
		if verbosity>0:
			"Constructing subhalo selection cut"

		mask = np.zeros(baryons.size)
		for name in cuts.keys():
			lower, upper = cuts[name].split()

			if verbosity>0:
				print name, lower, upper

			if '_dm' in name:
				sel = (dm[name.replace('_dm', '')]<float(upper)) & (dm[name.replace('_dm', '')]>float(lower))
				sel = sel & np.isfinite(dm[name.replace('_dm', '')])
			elif '_baryon' in name:
				sel = (baryons[name.replace('_baryon', '')]<float(upper)) & (baryons[name.replace('_baryon', '')]>float(lower))
				sel = sel & np.isfinite(baryons[name.replace('_baryon', '')])
			else:
				sel = (baryons[name]<float(upper)) & (baryons[name]>float(lower))
				sel = sel & np.isfinite(baryons[name])
			mask[sel] = 1

		self.mask = mask.astype(bool)

		if (verbosity>0):
			print 'Cuts leave %3.3f/%3.3f M galaxies (%3.3f percent)'%(len(self.mask[self.mask])*1./1e6, len(self.mask)*1./1e6, 100.*len(self.mask[self.mask])/len(self.mask))

		return self.mask

	def occupation_statistics(self):
		i,nocc = np.unique(self.array['halo_id'], return_counts=True)

		# There must be a better way of doing this.
		for (j,N)in zip(i,nocc):
			sel = (self.array['halo_id']==j)
			self.array['nocc'][sel] = N

		return

	def find_cs_flag(self):
		for ih in np.unique(self.array['halo_id']):

			rh = self.array['rh'][(self.array['halo_id']==ih)]
			Rmin = rh.min()
			ic = np.argwhere(rh==Rmin)
			self.array['central'][ic] = 1

	def export(self, outpath):

		# Check whether the file exists already, and if so remove the old version
		if os.path.exists(outpath):
			import pdb ; pdb.set_trace()
			os.system('rm %s'%outpath)

		print 'Saving file to %s'%outpath

		outfits = fi.FITS(outpath,'rw')
		outfits.write(self.array)
		outfits.close()

		print 'Done.'

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--verbosity', '-v', type=int, action='store', default=1)
parser.add_argument('--config', '-c', type=str, action='store')
args = parser.parse_args()

config = yaml.load(open(args.config))


# This is a table of subhalo data that Ananth compiled at some point
# Contains the basic positions (defined by the potential well), masses and particle numbers 
root_folder='/physics/yfeng1/mb2'
snapshot='0%d'%config['snapshot']

snap = SnapDir(snapshot, root_folder)
h = snap.readsubhalo()


# Read in the data
if args.verbosity>0:
	print 'Reading baryon information from ', config['baryon_shapes']
baryons = fi.FITS(config['baryon_shapes'])['baryons'].read()

if args.verbosity>0:
	print 'Reading dark matter information from ', config['dm_shapes']
dm = fi.FITS(config['dm_shapes'])['dm'].read()



# Now create an array in which to store the required columns
cat = catalogue()
# Impose the selection cuts, as specified in the configuration file
mask = cat.build_selection_mask(config['cuts'], baryons, dm, verbosity=args.verbosity)

cat.array = np.zeros(baryons[mask].size, dtype=dt)

# The ids, masses and particle numbers are fairly straightforwards to obtain 
i = np.arange(0, baryons.size, 1)
cat.array['object_id'] = i[mask]
cat.array['baryon_mass'] = h['massbytype'].T[1][mask]
cat.array['matter_mass'] = h['massbytype'].T[4][mask]
cat.array['npart_dm'] = h['lenbytype'].T[1][mask]
cat.array['npart_baryon'] = h['lenbytype'].T[4][mask]

# Also the positions
cat.array['x'] = h['pos'].T[0][mask]/1e3
cat.array['y'] = h['pos'].T[1][mask]/1e3
cat.array['z'] = h['pos'].T[2][mask]/1e3

# Copy over the 3D shape vectors
cat.array['a1'] = baryons['c1'][mask]
cat.array['a2'] = baryons['c2'][mask]
cat.array['a3'] = baryons['c3'][mask]
cat.array['b1'] = baryons['b1'][mask]
cat.array['b2'] = baryons['b2'][mask]
cat.array['b3'] = baryons['b3'][mask]
cat.array['c1'] = baryons['a1'][mask]
cat.array['c2'] = baryons['a2'][mask]
cat.array['c3'] = baryons['a3'][mask]

cat.array['a1_dm'] = dm['c1'][mask]
cat.array['a2_dm'] = dm['c2'][mask]
cat.array['a3_dm'] = dm['c3'][mask]
cat.array['b1_dm'] = dm['b1'][mask]
cat.array['b2_dm'] = dm['b2'][mask]
cat.array['b3_dm'] = dm['b3'][mask]
cat.array['c1_dm'] = dm['a1'][mask]
cat.array['c2_dm'] = dm['a2'][mask]
cat.array['c3_dm'] = dm['a3'][mask]

# Now calculate the projected ellipticities
phi = np.arctan2(cat.array['a2'], cat.array['a1'])
q = np.sqrt(baryons['lambda2'][mask]/baryons['lambda1'][mask])
e = (q-1)/(q+1)
e[np.invert(np.isfinite(e))] = 0.0
cat.array['e1'] = e * np.cos(2*phi)
cat.array['e2'] = e * np.sin(2*phi)


# And the same for the projected DM shapes
phi = np.arctan2(cat.array['a2_dm'], cat.array['a1_dm'])
q = np.sqrt(dm['lambda2'][mask]/dm['lambda1'][mask])
e = (q-1)/(q+1)
e[np.invert(np.isfinite(e))] = 0.0
cat.array['e1_dm'] = e * np.cos(2*phi)
cat.array['e2_dm'] = e * np.sin(2*phi)

# Now find the associated halo per galaxy
# We'll need to read out this information from the database
if config['include']['halo_matching']:
	halos = halo_wrapper(config['snapshot'], verbosity=args.verbosity)
	R, ind, info = halos.populate(cat.array)

	cat.array['halo_id'] = halos.ids[ind]
	cat.array['rh'] = R

	# Working out the halo occupation is slightly more fiddly
	cat.occupation_statistics()

# Now save the compiled subhalo data as a FITS file
cat.export(config['output'])

import pdb ; pdb.set_trace()
