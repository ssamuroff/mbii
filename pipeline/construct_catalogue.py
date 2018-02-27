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
      ('vrot', float), # Rotational velocity about the subhalo centre
      ('sigma', float), # Velocity dispersion  
      ('central', int), # Flag to indicate central galaxies
      ('rh', float), # Radial distance from the centre of the halo
      ('most_massive', int), # Binary flag, 1 if subhalo is has the highest DM mass in its group 
      ('spatial_central', float), # Bianry flag, 1 if subhalo is the closest to the potential minimum of the host halo 
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

		# Setup the database connection and query for the centroids of all the halos
		c = db.cursor()
		sql = 'SELECT x,y,z,groupId,mass,len FROM subfind_halos WHERE snapnum=85;'
		if verbosity>0:
			print 'Submitting query...'
		c.execute(sql)
		self.info = fromarrays(np.array(c.fetchall()).squeeze().T,names='x,y,z,groupId,mass,len')

		# Do the same for groups
		sql = 'SELECT x,y,z,groupId FROM subfind_groups WHERE snapnum=85;'
		if verbosity>0:
			print 'Submitting group query...'
		c.execute(sql)
		self.group_info = fromarrays(np.array(c.fetchall()).squeeze().T,names='x,y,z,groupId')

		# Convert to Mpc h^-1
		for name in ['x', 'y', 'z']:
			self.info[name]/=1e3
			self.group_info[name]/=1e3

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

		if self.verbosity>0:
			print 'Querying tree'

		# Query it for each 3D galaxy position
		xyz = np.array([data['x'], data['y'], data['z']])
		R,ind = tree.query(xyz.T, k=1)

		return R, ind, self.info[ind]

	def groups(self, data):
		# Match up each galaxy (subhalo) to the nearest group centroid
		if self.verbosity>0:
			print 'Building group KD tree'

		# Build the tree
		xyz0 = np.array([self.group_info['x'], self.group_info['y'], self.group_info['z']])
		tree = sps.KDTree(xyz0.T)

		if self.verbosity>0:
			print 'Querying tree'

		# Query it for each 3D galaxy position
		xyz = np.array([data['x'], data['y'], data['z']])
		R,ind = tree.query(xyz.T, k=1)

		return R, ind, self.group_info[ind]


class catalogue:
	def __init__(self):
		print 'Set up subhalo catalogue'

	def build_selection_mask(self, cuts, baryons, dm, h, verbosity=1):
		if verbosity>0:
			"Constructing subhalo selection cut"

		self.verbosity = verbosity

		mask = np.ones(baryons.size)
		for name in cuts.keys():
			lower, upper = cuts[name].split()

			if (verbosity>0):
				print '%s < %s < %s'%(lower, name, upper)

			if ('npart_dm' in name):
				n = h['lenbytype'].T[0]
				sel = (n>float(upper)) | (n<float(lower))
				sel = sel | np.invert(np.isfinite(n))
			elif ('npart_baryon' in name):
				n = h['lenbytype'].T[4]
				sel = (n>float(upper)) | (n<float(lower))
				sel = sel | np.invert(np.isfinite(n))
			elif (name in ['x', 'y', 'z']):
				lookup = {'x':0, 'y':1, 'z':2}
				p = h['pos'].T[lookup[name]]
				sel = (p>float(upper)) | (p<float(lower))
				sel = sel | np.invert(np.isfinite(p))
			else:
				sel = (baryons[name]>float(upper)) | (baryons[name]<float(lower))
				sel = sel | np.invert(np.isfinite(baryons[name]))
			mask[sel] = 0

		self.mask = mask.astype(bool)

		if (verbosity>0):
			print 'Cuts leave %3.3f/%3.3f M galaxies (%3.3f percent)'%(len(self.mask[self.mask])*1./1e6, len(self.mask)*1./1e6, 100.*len(self.mask[self.mask])/len(self.mask))

		return self.mask

	def occupation_statistics(self):
		i,edges,nocc = np.unique(self.array['halo_id'], return_counts=True, return_index=True)

		edges = np.hstack((edges, np.array([-1])))

		# There must be a better way of doing this.
		for j,N in enumerate(nocc):
			self.array['nocc'][edges[j]:edges[(j+1)]] = N

		self.array['nocc'][-1] = nocc[-1]

		return

	def tenneti_info(self, h, mask):
		centralflag = np.fromfile('/home/rmandelb.proj/ananth/centralflag_085',dtype=np.uint32)

		contam_mask = np.isfinite(h['pos'].T[1]) & (h['pos'].T[0]<100000) & (h['pos'].T[1]<100000) & (h['pos'].T[2]<100000)
		cflag = centralflag[mask[contam_mask]]

		self.array['central'] = cflag

		return

	def calculate_galaxy_offsets(self):

		# Array of halos for which centroid positions are needed
		halo_ids = np.unique(self.array['halo_id'])

		R = np.zeros(self.array.size, dtype=[('x',float), ('y',float), ('z',float)]) - 9999.

		for ih in halo_ids:
			# Query the database to find the 3D position of the potential minimum
			info = utils.find_centre(ih)
			# Convert to h^-1 Mpc
			x0 = info['x']/1000
			y0 = info['y']/1000
			z0 = info['z']/1000

			# Select galaxies in this halo
			select = (self.array['halo_id']==ih)

			# Work out the offset in position
			R['x'][select] = x0 - self.array['x'][select]
			R['y'][select] = y0 - self.array['y'][select]
			R['z'][select] = z0 - self.array['z'][select]

		return R

	def find_cs_flag(self, group_data, h, subhalo_data_dm, subhalo_data_bar):
		# This is very similar to the neighbour search carried out earlier
		# This should produce two flag tables, with the criterion for being a central galaxy
		# either:  
		# (a) being the nearest object with star particles to the potential minimum
		# (b) being the most massive subhalo in the group 
		# The masks are shaped like the full subhalo data, with mass cuts etc to be applied later

		print 'Constructing central flag table'

		groups = np.unique(group_data['groupId'])
		ngrp = len(groups)
		nbar = subhalo_data_bar['npart']
		cflag = np.zeros(subhalo_data_dm['x'].size)
		massflag = np.zeros(subhalo_data_dm['x'].size)
		i0=0

		for i, g in enumerate(group_data):
		    # Match up each galaxy (subhalo) to the nearest group centroid
		    print 'Building KD tree'

		    # Build the tree using galaxies in the group
		    # The cuts here reject
		    # (a) subhalos not associated with group i
		    # (b) subhalos with which no star particles are associated
		    # (c) subhalos with non-finite positions, or positions outside the simulation box
		    select = (h['groupid']==g['groupId']) & (nbar>0) & (h['pos'].T[0]/1000>0) & (h['pos'].T[0]/1000<100) & (h['pos'].T[1]/1000>0) & (h['pos'].T[1]/1000<100) & (h['pos'].T[2]/1000>0) & (h['pos'].T[2]/1000<100)

		    xyz0 = np.array([h['pos'].T[0][select]/1000, h['pos'].T[1][select]/1000, h['pos'].T[2][select]/1000])
		    tree = sps.KDTree(xyz0.T)

		    print 'Querying tree for group %d'%i

		    # Query it for the 3D halo centroid (one 3 vector)
		    xyz = np.array([g['x'], g['y'], g['z']])
		    R,ind = tree.query(xyz.T, k=1)

		    subflags = np.zeros(xyz0[0].size)
		    subflags[ind]+=1
		    cflag[select] = subflags

		    # Also flag the most massive subhalo in the group (by matter, not baryonic mass)
		    msubflags = np.zeros(xyz0[0].size)
		    msubflags[(h['lenbytype'].T[0][select]==h['lenbytype'].T[0][select].max())]+=1
		    massflag[select]=msubflags

		return cflag. massflag

		
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

# Do the central/satellite flagging
halos = halo_wrapper(config['snapshot'], verbosity=args.verbosity)
spatial_cflag, mass_cflag = cat.find_cs_flag(halos.group_info, h, dm, baryons)

# Impose the selection cuts, as specified in the configuration file
mask = cat.build_selection_mask(config['cuts'], baryons, dm, h, verbosity=args.verbosity)
mask = mask & np.isfinite(h['pos'].T[0]) & np.isfinite(h['pos'].T[1]) & np.isfinite(h['pos'].T[2])

cat.array = np.zeros(baryons[mask].size, dtype=dt)

cat.array['most_massive'] = mass_cflag[mask]
cat.array['spatial_central'] = spatial_cflag[mask]

# The ids, masses and particle numbers are fairly straightforwards to obtain 
i = np.arange(0, baryons.size, 1)
cat.array['object_id'] = i[mask]
cat.array['baryon_mass'] = h['massbytype'].T[4][mask]
cat.array['matter_mass'] = h['massbytype'].T[1][mask]
cat.array['npart_dm'] = h['lenbytype'].T[1][mask]
cat.array['npart_baryon'] = h['lenbytype'].T[4][mask]

# Also the positions
cat.array['x'] = h['pos'].T[0][mask]/1e3
cat.array['y'] = h['pos'].T[1][mask]/1e3
cat.array['z'] = h['pos'].T[2][mask]/1e3

cat.array['vrot'] = h['vcirc'][mask]
cat.array['sigma'] = h['vdisp'][mask]

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
phi = np.arctan2(baryons['a2_2d'][mask], baryons['a1_2d'][mask])
q = np.sqrt(baryons['lambda2_2d'][mask]/baryons['lambda1_2d'][mask])
e = (q-1)/(q+1)
e[np.invert(np.isfinite(e))] = 0.0
cat.array['e1'] = e * np.cos(2*phi)
cat.array['e2'] = e * np.sin(2*phi)


# And the same for the projected DM shapes
phi = np.arctan2(dm['a2_2d'][mask], dm['a1_2d'][mask])
q = np.sqrt(dm['lambda2_2d'][mask]/dm['lambda1_2d'][mask])
e = (q-1)/(q+1)
e[np.invert(np.isfinite(e))] = 0.0
cat.array['e1_dm'] = e * np.cos(2*phi)
cat.array['e2_dm'] = e * np.sin(2*phi)

# Now find the associated halo per galaxy
# We'll need to read out this information from the database
if config['include']['halo_matching']:
	halos = halo_wrapper(config['snapshot'], verbosity=args.verbosity)
	cat.find_cs_flag(halos.group_info)
	R, ind, info = halos.populate(cat.array)

	cat.array['halo_id'] = info['groupId']

	Rg, indg, infog = halos.groups(cat.array)
	cat.array['rh'] = Rg



else:
	ref = fi.FITS('/physics2/ssamurof/massive_black_ii/cats/base_subhalo_shapes-v2.fits')[1].read()
	# Otherwise, just copy over the relevant columns from the base catalogue
	for name in ['central', 'nocc', 'rh', 'halo_id']:
		cat.array[name] = ref[name]

# Working out the halo occupation is slightly more fiddly
cat.occupation_statistics()
cat.tenneti_info(h, mask)
	
#cat.calculate_galaxy_offsets()

import pdb ; pdb.set_trace
# Now save the compiled subhalo data as a FITS file
cat.export(config['output'])

import pdb ; pdb.set_trace()

