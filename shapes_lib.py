import numpy as np
import copy
import fitsio as fi
import pylab as plt
plt.switch_backend('agg')
import numpy as np
import mbii.lego_tools as utils
from numpy.core.records import fromarrays

import numpy as np
from readsubhalo import *
from mbii.readsubhalo import *
import fitsio as fi

colnames = {'SubhaloMassType':'massbytype', 'SubhaloLenType':'lenbytype', 'SubhaloGrNr':'groupid', 'pos':'SubhaloPos'} 

def build_matrix(pos, reduced=False):
	"""Construct a 3x3 inertia tensor of particles in a subhalo.
	Requires an array of particle positions relative to the subhalo centroid."""

	# Normalisation factor
	norm = 1.0 * len(pos[0])

	#Select the appropriate per-particle weights
	if reduced:
		wt = np.sum(pos*pos, axis=0)
		select = (wt!=0)
	else:
		wt = np.ones(pos[0].size)
		select = wt.astype(bool)

	# Put together the 3x3 matrix
	tensor = np.zeros((3,3))
	tensor[0, 0] = np.dot( pos[0][select], pos[0][select]/wt[select])/norm
	tensor[1, 1] = np.dot( pos[1][select], pos[1][select]/wt[select])/norm
	tensor[2, 2] = np.dot( pos[2][select], pos[2][select]/wt[select])/norm
	tensor[1, 0] = np.dot( pos[1][select], pos[0][select]/wt[select])/norm
	tensor[2, 0] = np.dot( pos[2][select], pos[0][select]/wt[select])/norm
	tensor[2, 1] = np.dot( pos[2][select], pos[1][select]/wt[select])/norm
	tensor[0, 2] = tensor[2, 0]
	tensor[0, 1] = tensor[1, 0]
	tensor[1, 2] = tensor[2, 1]

	return tensor

def compute(i, k, x, subhalo_centroids, eigvalues, eigvalues_2d, eigvectors, eigvectors_2d, length, reduced=False):
	norm = 1.0 * len(x.T[0])

	# Check each coordinate, as calculated, is finite
	# If not, use the weighted mean position of the dark matter particles instead
	for j in [0,1,2]:
		if not np.isfinite(subhalo_centroids[j][i]):
			subhalo_centroids[j][i] = np.sum(x.T[j])/norm

	# Use the centre of the subhalo's potential well as the centroid here
	# Work out the radial position of each particle from the centroid
	pos = np.array([ x.T[0] - subhalo_centroids[0][i], 
		             x.T[1] - subhalo_centroids[1][i], 
		             x.T[2] - subhalo_centroids[2][i] ])

	tens = build_matrix(pos, reduced=reduced)

	# Compute the eigenvalues of the halos and store the outputs
	w, v = np.linalg.eigh(tens)
	w2d, v2d = np.linalg.eigh(tens[:2,:2])

	eigvalues[k,i] = w
	eigvalues_2d[k,i] = w2d
	eigvectors[k,i] = v
	eigvectors_2d[k,i] = v2d
	length[k,i] = norm

	return eigvalues, eigvalues_2d, eigvectors, eigvectors_2d, length

def export(inclusion_threshold, eigvalues, eigvalues_2d, eigvectors_2d, eigvectors, subhalo_centroids, length, reduced=False, snapshot=85, simulation='massiveblackii', dirname='', rank=0):
	"""Sorry."""

	# Decide on the filename
 	if reduced:
 		filename = '%s/%s-subhalo_cat_reduced-nthreshold%d-snapshot%d-proj+3d-%d.fits'%(dirname, simulation, inclusion_threshold, snapshot, rank)
 	else:
 		filename = '%s/%s-subhalo_cat-nthreshold%d-snapshot%d-proj+3d-%d.fits'%(dirname, simulation, inclusion_threshold, snapshot, rank)

 	print "Saving output to %s"%filename
 	out = fi.FITS(filename,'rw')

 	dt = [('x', float), ('y', float), ('z', float), ('npart',float), ('lambda1', float), ('lambda2', float), ('lambda3', float), ('a1', float), ('a2', float), ('a3', float), ('b1', float), ('b2', float), ('b3', float), ('c1', float), ('c2', float), ('c3', float), ('lambda1_2d', float), ('lambda2_2d', float), ('a1_2d', float), ('a2_2d', float), ('b1_2d', float), ('b2_2d', float)]
 	dat=np.zeros(eigvalues[0].T[0].size, dtype=dt)

 	dat['lambda1'] = eigvalues[0].T[0]
 	dat['lambda2'] = eigvalues[0].T[1]
 	dat['lambda3'] = eigvalues[0].T[2]
 	dat['lambda1_2d'] = eigvalues_2d[0].T[0]
 	dat['lambda2_2d'] = eigvalues_2d[0].T[1]
 	dat['a1_2d'] = eigvectors_2d[0].T[0,0]
 	dat['a2_2d'] = eigvectors_2d[0].T[0,1]
 	dat['b1_2d'] = eigvectors_2d[0].T[1,0]
 	dat['b2_2d'] = eigvectors_2d[0].T[1,1]
 	dat['a1'] = eigvectors[0].T[0,0]
 	dat['a2'] = eigvectors[0].T[0,1]
 	dat['a3'] = eigvectors[0].T[0,2]
 	dat['b1'] = eigvectors[0].T[1,0]
 	dat['b2'] = eigvectors[0].T[1,1]
 	dat['b3'] = eigvectors[0].T[1,2]
 	dat['c1'] = eigvectors[0].T[2,0]
 	dat['c2'] = eigvectors[0].T[2,1]
 	dat['c3'] = eigvectors[0].T[2,2]
 	dat['x'] = subhalo_centroids[0]
 	dat['y'] = subhalo_centroids[1]
 	dat['z'] = subhalo_centroids[2]
 	dat['npart'] = length[0]

 	out.write(dat)
 	out[-1].write_key('EXTNAME', 'dm')

 	dt = [('x', float), ('y', float), ('z', float), ('npart',float), ('lambda1', float), ('lambda2', float), ('lambda3', float), ('a1', float), ('a2', float), ('a3', float), ('b1', float), ('b2', float), ('b3', float), ('c1', float), ('c2', float), ('c3', float), ('lambda1_2d', float), ('lambda2_2d', float), ('a1_2d', float), ('a2_2d', float), ('b1_2d', float), ('b2_2d', float)]
 	dat2=np.zeros(eigvalues[1].T[0].size, dtype=dt)

 	dat2['lambda1'] = eigvalues[1].T[0]
 	dat2['lambda2'] = eigvalues[1].T[1]
 	dat2['lambda3'] = eigvalues[1].T[2]
 	dat2['lambda1_2d'] = eigvalues_2d[1].T[0]
 	dat2['lambda2_2d'] = eigvalues_2d[1].T[1]
 	dat2['a1_2d'] = eigvectors_2d[1].T[0,0]
 	dat2['a2_2d'] = eigvectors_2d[1].T[0,1]
 	dat2['b1_2d'] = eigvectors_2d[1].T[1,0]
 	dat2['b2_2d'] = eigvectors_2d[1].T[1,1]
 	dat2['a1'] = eigvectors[1].T[0,0]
 	dat2['a2'] = eigvectors[1].T[0,1]
 	dat2['a3'] = eigvectors[1].T[0,2]
 	dat2['b1'] = eigvectors[1].T[1,0]
 	dat2['b2'] = eigvectors[1].T[1,1]
 	dat2['b3'] = eigvectors[1].T[1,2]
 	dat2['c1'] = eigvectors[1].T[2,0]
 	dat2['c2'] = eigvectors[1].T[2,1]
 	dat2['c3'] = eigvectors[1].T[2,2]
 	dat2['x'] = subhalo_centroids[0]
 	dat2['y'] = subhalo_centroids[1]
 	dat2['z'] = subhalo_centroids[2]
 	dat2['npart'] = length[1]

 	out.write(dat2)
 	out[-1].write_key('EXTNAME', 'baryons')
 	out.close()

def load_single_subhalo(i, simulation, snapshot, root):
	if (simulation.lower()=='massiveblackii'):

		snap = SnapDir('0%d'%snapshot, root)
		h = snap.readsubhalo()

		# Load the positions and masses of the constituent particles
		#print 'Loading dark matter particles'
		x = snap.load(1, 'pos', h)[i]

		#print 'Loading star particles'
		xb = snap.load(4, 'pos', h)[i]

	elif (simulation.lower()=='illustris'):
		import illustris_python as il

		#print 'Loading dark matter particles'
		x = il.snapshot.loadSubhalo(root, snapshot, i ,'dm',['Coordinates'])
		#print 'Loading star particles'
		xb = il.snapshot.loadSubhalo(root, snapshot, i ,'stars',['Coordinates'])

	else:
		raise ValueError('Unknown simulation. Sorry.')

	return x, xb


def read_subhalo_data(simulation, snapshot, root, columns=['SubhaloPos']):
	"""Read in the particle information from either MBII or Illustris.
	   The formats of the data on disc are slightly different, and so
	   some fiddling is needed to get the catalogues into something
	   the shape code can use equivalently."""

	if (simulation.lower()=='massiveblackii'):

		snap = SnapDir('0%d'%snapshot, root)
		h = snap.readsubhalo()

		# Load the positions and masses of the constituent particles
		print 'Loading dark matter particles'
		x = snap.load(1, 'pos', h)

		print 'Loading star particles'
		xb = snap.load(4, 'pos', h)

		subhalo_positions = h['pos']

	elif (simulation.lower()=='illustris'):
		import illustris_python as il
		x = []
		m = []
		xb = []
		mb = []

		subhalo_positions = il.groupcat.loadSubhalos(root, snapshot, fields=columns)

	else:
		raise ValueError('Unknown simulation. Sorry.')

	return subhalo_positions, x, xb

def read_subhalo_data_all(simulation, snapshot, root):

	if (simulation.lower()=='massiveblackii'):

		snap = SnapDir('0%d'%snapshot, root)
		h = snap.readsubhalo()

		ms = h['massbytype'].T[4]
		md = h['massbytype'].T[1]
		mg = h['massbytype'].T[0]
		nd = h['lenbytype'].T[1]
		ns = h['lenbytype'].T[4]
		vdisp = h['vdisp'] 
		spin = h['vcirc']
		ig = h['groupid']


		pos = h['pos']

	elif (simulation.lower()=='illustris'):
		import illustris_python as il

		info = il.groupcat.loadSubhalos(root, snapshot)

		pos = info['SubhaloPos']

		ms = info['SubhaloMassType'].T[4]
		md = info['SubhaloMassType'].T[1]
		mg = info['SubhaloMassType'].T[0]
		nd = info['SubhaloLenType'].T[1]
		ns = info['SubhaloLenType'].T[4]

		vdisp = info['SubhaloVelDisp']
		spin = info['SubhaloSpin']
		ig = info['SubhaloGrNr']

	else:
		raise ValueError('Unknown simulation. Sorry.')

	return ig, pos, md, ms, mg, nd, ns, vdisp, spin

def compute_inertia_tensors(options, rank, size, reduced=False, inclusion_threshold=5, snapshot=85, savedir=''):
 	""" Compute the intertia tensors for all subhalos. Do the calculation twice, for the
 	the dark matter and stellar component. Flatten the results and save as columns in a FITS file."""
 	print 'Using reduced (distance weighted) tensors', 'yes'*int(reduced), 'no'*int(np.invert(reduced) )

	# This is horrible. Sorry. 
	# Ultimately I'd like to eliminate the last bits of code inherited from Ananth,
	# but digging into how the SubFind outputs are accessed requires a non-trivial amount of work.  

	# Read the subhalo information

	simulation = options['simulation']
	snapshot = options['catalogues']['snapshot']
	root_folder = options['root_folder']
	h, x, xb = read_subhalo_data(simulation, snapshot, root_folder)
	
    
	eigvectors = np.zeros((2, len(h), 3, 3))
	eigvalues  = np.zeros((2, len(h), 3))
	eigvectors_2d = np.zeros((2, len(h), 2, 2))
	eigvalues_2d  = np.zeros((2, len(h), 2))
	length  = np.zeros((2, len(h)))

	subhalo_centroids = np.array(h.T)
 
	# Will compute an inertia tensor per subhalo
	for i in range(len(h)):

		# MPI handline
		if i%size!=rank:
			continue

		if (i%100 ==0):
			print "Done %d samples"%i

		if (len(x)>0):
			particle_positions_dm = x[i]
			particle_positions_b = xb[i]
		else:
			particle_positions_dm, particle_positions_b = load_single_subhalo(i, simulation, snapshot, root_folder)

		# Reject subhalos with less than some threshold occupation number
		if (len(particle_positions_dm) < inclusion_threshold):
			pass
		else: 
			eigvalues, eigvalues_2d, eigvectors, eigvectors_2d, length = compute(i, 0, particle_positions_dm, subhalo_centroids, eigvalues, eigvalues_2d, eigvectors, eigvectors_2d, length, reduced=reduced)
			
		if (len(particle_positions_b) < inclusion_threshold):
			pass
		else:
			eigvalues, eigvalues_2d, eigvectors, eigvectors_2d, length = compute(i, 1, particle_positions_b, subhalo_centroids, eigvalues, eigvalues_2d, eigvectors, eigvectors_2d, length, reduced=reduced)

	export(inclusion_threshold, eigvalues, eigvalues_2d, eigvectors_2d, eigvectors, subhalo_centroids, length, reduced=reduced, snapshot=snapshot, simulation=simulation, dirname=savedir, rank=rank)
	print 'Done'

	return
 	




def compute_spin(options, rank, size, inclusion_threshold=1, component='baryons', nsubhalo=-1):
    """ Compute the angular momentum for all subhalos for
    the dark matter and stellar components and saves the output
    as a FITS file.
    """

    indices = {'baryons':4, 'dm':1}

    # Read the subhalo information
    h = snap.readsubhalo()
    # Load the positions and masses of the constituent particles
    x = snap.load(indices[component], 'pos', h)
    m = snap.load(indices[component], 'mass', h)
    v = snap.load(indices[component], 'vel', h)
    
    Vr = np.zeros(len(h))
    Sigma = np.zeros(len(h))
    spin  = np.zeros((len(h), 3))

 
    # Will compute for each halo the inertia tensor
    if nsubhalo<0:
        nsubhalo=len(h)
        print 'Will process all (%d) subhalos'%nsubhalo

    print 'Inclusion Threshold : %d particles'%inclusion_threshold

    for i in range(nsubhalo):
        if i%100 ==0:
            print "Done %d samples"%i

        if len(x[i]) < inclusion_threshold:
            pass # print "Halo %d is empty of dark matter"%i
        else:
            weights = np.ones(len(x[i]))
            normFactor = np.double(len(x[i].T[0]))
            x0 = x[i].T[0] - np.dot(x[i].T[0],weights)/normFactor
            x1 = x[i].T[1] - np.dot(x[i].T[1],weights)/normFactor
            x2 = x[i].T[2] - np.dot(x[i].T[2],weights)/normFactor

            r = np.array([x0,x1,x2])
            # Calaculate the angular momentum vectors for individual particles 
            l = [M*np.cross(X,V) for (M,X,V) in zip(m[i],r.T,v[i])]
            # Sum them to get a three vector halo spin
            L = np.sum(l,axis=0)

            # Rotate the Cartesian velocity into a frame defined by the angular momentum vector
            phi_xy = np.arctan2(L[1],L[0])
            phi_yz = np.arctan2(L[1],L[2])
            phi_xz = np.arctan2(L[0],L[2])

            Rx = np.array([[1,0,0], [0,np.cos(phi_yz), -np.sin(phi_yz)], [0,np.sin(phi_yz),np.cos(phi_yz)]])
            Ry = np.array([[np.cos(phi_xz),0,np.sin(phi_xz)], [0,1,0], [-np.sin(phi_xz),0,np.cos(phi_xz)]])
            Rz = np.array([[np.cos(phi_xy), -np.sin(phi_xy), 0], [np.sin(phi_xy),np.cos(phi_xy),0], [0,0,1]])

            # Apply the rotations about each axis in turn
            # The result is a three component vector vx,vy,vz with vz aligned with the spin axis
            vec = np.dot(Rx,v[i].T)
            vec = np.dot(Ry,vec)
            vec = np.dot(Rz,vec)

            # This is the same for the particle position vectors
            rvec = np.dot(Rx,r)
            rvec = np.dot(Ry,rvec)
            rvec = np.dot(Rz,rvec)

            # Work out the 2D position and velocity magnitudes in the plane of the disc
            #r0 = np.sqrt(rvec[0]*rvec[0] + rvec[1]*rvec[1])
            t1 = np.arctan2(rvec[1],rvec[0])
            #v0 = np.sqrt(vec[0]*vec[0] + vec[1]*vec[1])
            t2 = np.arctan2(vec[1],vec[0])
            theta_rv = t2-t1

            # The angle between the separation vector and the velocity vector
            #theta_rv = (rvec[0]*vec[0] + rvec[1]*vec[1]) /r0 /v0
            #theta_rv = np.arccos(theta_rv)

            vrot = v0 * np.sin(theta_rv)

            #import pdb ; pdb.set_trace()

            sigma = np.array([v[i].T[0].std(), v[i].T[1].std(), v[i].T[2].std() ])
            sigma = np.sqrt(sum(sigma*sigma))

            #theta = np.arctan2(rvec[1],rvec[0])
            #vrot = vec[1]*np.cos(theta) - vec[0]*np.sin(theta)
            v_sigma = np.mean(vrot)/sigma

            #import pdb ; pdb.set_trace()

            spin[i] = L
            Vr[i] = np.mean(vrot)
            Sigma[i] = sigma

    print "Saving output"
    #np.save("eigenvaluesb.npy", eigvalues)
    #np.save("eigenvectorsb.npy", eigvectors)
    #np.save("centroids.npy", centroids)

    out = fi.FITS('/home/ssamurof/massive_black_ii/subhalo_spin-%s-nthreshold%d-nsubhalo%d.fits'%(component, inclusion_threshold, nsubhalo),'rw')
    

    dat=np.zeros(len(h), dtype=[ ('vrot', float), ('sigma', float), ('spinx', float), ('spiny', float), ('spinz', float)])

    dat['vrot'] = Vr
    dat['sigma'] = Sigma
    dat['spinx'] = spin.T[0]
    dat['spiny'] = spin.T[1]
    dat['spinz'] = spin.T[2]

    out.write(dat)
    out[-1].write_key('EXTNAME', component)

    out.close()


