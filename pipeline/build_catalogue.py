import numpy as np
import argparse
import yaml
import fitsio as fi
import scipy.spatial as sps
from mbii.readsubhalo import *
from mbii.properties import *
import mbii.shapes_lib as slib
import pymysql as mdb
import numpy as np
import mbii.symmetrise_lib as lib
from numpy.core.records import fromarrays


dt = [
      ('object_id', int), # unique ID for each galaxy 
      ('gas_mass', float), # What it sounds like, in units of 10^10 h^-1 M_*
      ('matter_mass', float),  # What it sounds like, in units of 10^10 h^-1 M_*
      ('stellar_mass', float),  # What it sounds like, in units of 10^10 h^-1 M_*
      ('npart_dm', int), # Number of DM particles in the subhalo
      ('npart_baryon', int), # Number of star particles in the subhalo
      ('e1', float), ('e2', float), # Projected ellipticities
      ('a1', float), ('a2', float), ('a3', float), # Normalised 3 vector along the major axis 
      ('b1', float), ('b2', float), ('b3', float),
      ('c1', float), ('c2', float), ('c3', float), # Normalised 3 vector along the major axis 
      ('lambda_a', float), ('lambda_b', float), ('lambda_c', float), # Eigenvalues of the 3D inertia tensor
      ('lambda_a_proj', float), ('lambda_b_proj', float), # Eigenvalues of the 2D inertia tensor 
      ('e1_dm', float), ('e2_dm', float), # Projected ellipticities, defined by the dark matter distribution
      ('a1_dm', float), ('a2_dm', float), ('a3_dm', float), # Normalised 3 vector along the major axis, defined by the dark matter distribution 
      ('b1_dm', float), ('b2_dm', float), ('b3_dm', float),
      ('c1_dm', float), ('c2_dm', float), ('c3_dm', float), # Normalised 3 vector along the major axis , defined by the dark matter distribution
      ('lambda_a_dm', float), ('lambda_b_dm', float), ('lambda_c_dm', float), # Eigenvalues of the 3D matter inertia tensor
      ('lambda_a_proj_dm', float), ('lambda_b_proj_dm', float), # Eigenvalues of the 2D matter inertia tensor 
      ('x', float), ('y', float), ('z', float), # 3D position in Cartesian coordinates, in Mpc h^-1
      ('xcent', float), ('ycent', float), ('zcent', float), # 3D position of the central galaxy in Cartesian coordinates, in Mpc h^-1
      ('vrot', float), # Rotational velocity about the subhalo centre
      ('sigma', float), # Velocity dispersion  
      #('central', int), # Flag to indicate central galaxies
      ('rh', float), # Radial distance from the centre of the halo
      ('most_massive', int), # Binary flag, 1 if subhalo is has the highest DM mass in its group 
      ('spatial_central', float), # Binary flag, 1 if subhalo is the closest to the potential minimum of the host halo
      ('hybrid_central', int), # Binary flag, 1 if subhalo is both in the 5 closest to the potential minimum and the 5 most massive objects in its host halo
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
		sql = 'SELECT x,y,z,groupId,mass,len,subfindId FROM subfind_halos WHERE snapnum=%d;'%snapshot
		if verbosity>0:
			print 'Submitting query...'
		c.execute(sql)
		self.info = fromarrays(np.array(c.fetchall()).squeeze().T,names='x,y,z,groupId,mass,len,subfindId')

		# Do the same for groups
		sql = 'SELECT x,y,z,groupId,subfindId FROM subfind_groups WHERE snapnum=%d;'%snapshot
		if verbosity>0:
			print 'Submitting group query...'
		c.execute(sql)
		self.group_info = fromarrays(np.array(c.fetchall()).squeeze().T,names='x,y,z,groupId,subfindId')

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

		# Build the tree of groups
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

	def build_selection_mask(self, cuts, q, verbosity=1):
		if verbosity>0:
			"Constructing subhalo selection cut"

		self.verbosity = verbosity

		mask = np.ones(baryons.size)
		for name in cuts.keys():
			lower, upper = cuts[name].split()

			if (verbosity>0):
				print '%s < %s < %s'%(lower, name, upper)

			sel = (q[name]>float(upper)) | (q[name]<float(lower))
			sel = sel | np.invert(np.isfinite(q[name]))
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
		# Not using these any more, because Ananth's central flags are crap
		# where
		# crap == wrong
		# (the flagging appears to have been done after mass cuts have been applied)
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

		# First do the KD tree calculation on the global level (ie not split by group)
		# This should hopefully be slightly faster than repeating the calculation for every halo
		# though maybe not...
		select = (nbar>0) & (h['pos'].T[0]/1000>0) & (h['pos'].T[0]/1000<100) & (h['pos'].T[1]/1000>0) & (h['pos'].T[1]/1000<100) & (h['pos'].T[2]/1000>0) & (h['pos'].T[2]/1000<100)

		xyz0 = np.array([h['pos'].T[0][select]/1000, h['pos'].T[1][select]/1000, h['pos'].T[2][select]/1000])

		print '---- Building KD tree from %d objects'%h['groupid'][select].size
		tree = sps.KDTree(xyz0.T)

		print '---- Querying tree for all groups (%d)'%ngrp
		# Query it for the 3D halo centroid (one 3 vector per halo)
		xyz = np.array([group_data['x'], group_data['y'], group_data['z']])
		R,ind = tree.query(xyz.T, k=1)

		# Remove any cases where the nearest galaxy is accociated with a different halo
		# Without this we'll be misclassifying a fair number of satellites as the centrals 
		# of neighbouring empty halos
		empty_mask = (h['groupid'][select][ind]==group_data['groupId'])

		# Store the results in the same shape as the subhalo table
		sub = np.zeros(cflag[select].shape)
		sub[ind[empty_mask]] = 1
		cflag[select] = sub

		# Do the same thing, but the other way around
		# to get a centroid-galaxy distance for each catalogue object
		print '---- Calculating galaxy offsets'
		inverted_tree = sps.KDTree(xyz.T)
		R_pergal,ind_pergal = inverted_tree.query(xyz0.T, k=1)

		Rp = np.zeros(subhalo_data_dm['x'].size)-1
		Rp[select] = R_pergal

		print 'Building mass flags.'
		Mm = np.array(h['lenbytype'].T[1]) # DM masses
		Mg = np.array(h['lenbytype'].T[0]) # gas masses
		gids = np.array(h['groupid'])

		# Attempting to do this without a for loop here
		# We'll need to assume the subhalos are sorted by mass within each group
		# which I _think_ is true for the historic data products
		ids,indices=np.unique(gids, return_index=True)
		massflag[indices]=1

#		for i, g in enumerate(group_data):
#		    # Match up each galaxy (subhalo) to the nearest group centroid
#
#		    # Build the tree using galaxies in the group
#		    # The cuts here reject
#		    # (a) subhalos not associated with group i
#		    # (b) subhalos with which no star particles are associated
#		    # (c) subhalos with non-finite positions, or positions outside the simulation box
#
		    #select = (h['groupid']==g['groupId']) & (nbar>0) & (h['pos'].T[0]/1000>0) & (h['pos'].T[0]/1000<100) & (h['pos'].T[1]/1000>0) & (h['pos'].T[1]/1000<100) & (h['pos'].T[2]/1000>0) & (h['pos'].T[2]/1000<100)
#
#		    # Also flag the most massive subhalo in the group (by matter, not baryonic mass)
#		    # Unfortunately this does require iteration over individual halos
#		    msubflags = np.zeros(h['groupid'][select].size)
#		    msubflags[(h['lenbytype'].T[0][select]==h['lenbytype'].T[0][select].max())]+=1
#		    massflag[select]=msubflags

		return cflag, massflag, Rp

		
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
root_folder=config['root_folder']
snapshot='0%d'%config['catalogues']['snapshot']


ig, pos, md, ms, mg, nd, ns, vspin, vdisp = slib.read_subhalo_data_all(config['simulation'], int(snapshot), root_folder)


#snap = SnapDir(snapshot, root_folder)
#h = snap.readsubhalo()


# Read in the data
if (config['catalogues']['shapes_method']=='reduced_inertia_tensor'):
	shapes_filename = '%s/%s-subhalo_cat_reduced-nthreshold%d-snapshot%d-proj+3d.fits'%( config['catalogues']['shapes_dir'], config['simulation'], 5, int(snapshot))
elif (config['catalogues']['shapes_method']=='inertia_tensor'):
	shapes_filename = '%s/%s-subhalo_cat-nthreshold%d-snapshot%d-proj+3d.fits'%( config['catalogues']['shapes_dir'], config['simulation'], 5, int(snapshot))
else:
	raise ValueError("Unknown shape method. Please choose from 'reduced_inertia_tensor' and 'inertia_tensor'.")

if args.verbosity>0:
	print 'Reading subhalo information from ', shapes_filename

dm = fi.FITS(shapes_filename)['dm'].read()
baryons = fi.FITS(shapes_filename)['baryons'].read()

# Now create an array in which to store the required columns
cat = catalogue()

# Impose the selection cuts, as specified in the configuration file
h = np.zeros(md.size, dtype=[('npart_dm',int), ('npart_baryon', int), ('x', float), ('y', float), ('z', float), ('groupid', int) ])
h['x'] = pos.T[0]
h['y'] = pos.T[1]
h['z'] = pos.T[2]
h['npart_baryon']=ns
h['npart_dm']=nd
mask = cat.build_selection_mask(config['catalogues']['cuts'], h, verbosity=args.verbosity)
mask = mask & np.isfinite(pos.T[0]) & np.isfinite(pos.T[1]) & np.isfinite(pos.T[2])

cat.array = np.zeros(baryons[mask].size, dtype=dt)

print 'Building mass flags.'
massflag = np.zeros(dm['x'].size)
Mm = np.array(md) # DM masses
gids = np.array(ig)
ids,indices=np.unique(gids, return_index=True)
massflag[indices]=1

cat.array['most_massive'] = massflag[mask]

# The ids, masses and particle numbers are fairly straightforwards to obtain 
i = np.arange(0, baryons.size, 1)
cat.array['object_id'] = i[mask]
cat.array['stellar_mass'] = ms[mask]
cat.array['matter_mass'] = md[mask]
cat.array['gas_mass'] = mg[mask]
cat.array['npart_dm'] = nd[mask]
cat.array['npart_baryon'] = ns[mask]

# Also the positions
cat.array['x'] = pos.T[0][mask]/1e3
cat.array['y'] = pos.T[1][mask]/1e3
cat.array['z'] = pos.T[2][mask]/1e3

cat.array['vrot'] = vspin[mask]
#cat.array['sigma'] = np.sqrt(np.sum([ vdisp.T[i0]**2 for i0 in [0,1,2] ], axis=0))

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


# Also copy over the information about the relative axis lengths
# of the visible and dark matter parts of each galaxy  
for (abc, num) in zip(['a','b','c'], [3,2,1]):
	cat.array['lambda_%c'%abc] = baryons['lambda%d'%num][mask]
	cat.array['lambda_%c_dm'%abc] = dm['lambda%d'%num][mask]

for (ab, num) in zip(['a','b'], [1,2]):
	cat.array['lambda_%c_proj'%ab] = baryons['lambda%d_2d'%num][mask]
	cat.array['lambda_%c_proj_dm'%ab] = dm['lambda%d_2d'%num][mask]


# Now calculate the projected ellipticities
for j in xrange(cat.array.size):
	cat.array = lib.project_ellipticities(j, cat.array, suffix='')

# And the same for the projected DM shapes
for j in xrange(cat.array.size):
	cat.array = lib.project_ellipticities(j, cat.array, suffix='dm')

# Now find the associated halo per galaxy
# We'll need to read out this information from the database
if config['catalogues']['halo_matching']:
	if (config['simulation']=='massiveblackii'):
		R, ind, info = halos.populate(cat.array)
		cat.array['halo_id'] = info['groupId']
	else:
		import illustris_python as il
		root='/nfs/nas-0-1/vat/Illustris-1'
		#sub = il.groupcat.loadSubhalos(root, config['catalogues']['snapshot'])
		cat.array['halo_id'] = ig[mask]
		h['groupid'] = ig

else:
	ref = fi.FITS('/physics2/ssamurof/massive_black_ii/cats/base_subhalo_shapes-v2.fits')[1].read()
	# Otherwise, just copy over the relevant columns from the base catalogue
	for name in ['central', 'nocc', 'rh', 'halo_id']:
		cat.array[name] = ref[name]

#import pdb ; pdb.set_trace()
#centres = [[lib.find_centre(g, snapshot=config['catalogues']['snapshot'])['x'], lib.find_centre(g, snapshot=config['catalogues']['snapshot'])['y'], lib.find_centre(g, snapshot=config['catalogues']['snapshot'])['z']] for g in np.unique(h['groupid']) ]

print 'Obtaining galaxy offsets'
for i, g in enumerate(cat.array['halo_id']):
	centroid = lib.find_centre(g, snapshot=config['catalogues']['snapshot'], simulation=config['simulation'])
	x0 = centroid['x']
	y0 = centroid['y']
	z0 = centroid['z']
	x = cat.array['x'][i]
	y = cat.array['y'][i]
	z = cat.array['z'][i]
	rh = np.sqrt((x0-x)**2 + (y0-y)**2 + (z0-z)**2)
	cat.array['rh'][i] = rh

print 'Building spatial central flag'
scflag=np.zeros(pos.T[0].size)
xcent = np.zeros(pos.T[0].size)
ycent = np.zeros(pos.T[1].size)
zcent = np.zeros(pos.T[2].size) 

ng=np.unique(h['groupid']).size
for i, g in enumerate(np.unique(ig)):
    if not (g in np.unique(cat.array['halo_id'])): continue 
    centroid = lib.find_centre(g, config['catalogues']['snapshot'], simulation=config['simulation'])
    x0 = centroid['x']
    y0 = centroid['y']
    z0 = centroid['z']
    hmask = (h['groupid']==g)
    x = pos.T[0][hmask]/1000
    y = pos.T[1][hmask]/1000
    z = pos.T[2][hmask]/1000
    rh = np.sqrt((x0-x)**2 + (y0-y)**2 + (z0-z)**2)
    flags = np.zeros(rh.size)
    fmask = np.isfinite(rh)
    flags[(rh==rh[fmask].min())] = 1
    scflag[hmask] = flags
    xcent[hmask] = np.atleast_1d(x[(rh==rh[fmask].min())])[0]
    ycent[hmask] = np.atleast_1d(y[(rh==rh[fmask].min())])[0]
    zcent[hmask] = np.atleast_1d(z[(rh==rh[fmask].min())])[0]
    print i, g 


cat.array['spatial_central'] = scflag[mask]
cat.array['xcent'] = xcent[mask]
cat.array['ycent'] = ycent[mask]
cat.array['zcent'] = zcent[mask]

# Working out the halo occupation is slightly more fiddly
cat.occupation_statistics()
#cat.tenneti_info(h, mask)
	
#cat.calculate_galaxy_offsets()

# Now save the compiled subhalo data as a FITS file
cat.export(config['catalogues']['postprocessed'])



