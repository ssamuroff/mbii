import numpy as np
import copy
import fitsio as fi
import pylab as plt
plt.switch_backend('agg')
import numpy as np
import mbii.lego_tools as utils
from numpy.core.records import fromarrays


def symmetrise_catalogue3(data=None, seed=4000, mask=None,filename='/home/ssamurof/massive_black_ii/subhalo_cat-nthreshold5.fits', savedir=None, replace=False, verbose=True, nmin=0, rank=0, size=1):
    """Takes a catalogue of subhalos with orientations, positions and host halo identifiers.
	   Cycles through the objects and applies a random rotation about the halo centre to each,
	   recalculating the position and shape vectors such that the relative orientation to
	   the centre is unchanged.
	   Then projects the 3D orientations along one axis of the simulation box and recomputes
	   the ellipticities.
	   Writes the result to a new FITS catalogue, identical to the input save for the shapes
	   and positions.
	"""

	# Set the random seed once at the top of the loop
    np.random.seed(seed)

    # Null mask if none is specified
    if mask is None:
    	mask = np.ones(data.size).astype(bool)

    # Safety catch
    print 'Will write to %s'%savedir
    if (savedir==filename):
    	print "WARNING: target and input file paths are the same"
    	import pdb ; pdb.set_trace()

    # Pull out the identifiers associating galaxies with halos
    group_ids = data['halo_id'][mask]
    groups = np.unique(group_ids)

    # Copy over the old data
    outdat = np.zeros(len(data), dtype=data.dtype)
    for name in outdat.dtype.names:
        outdat[name] = data[name]

    outdat = utils.add_col(outdat, 'x0', np.array([0]*outdat.size), dtype=float)
    outdat = utils.add_col(outdat, 'y0', np.array([0]*outdat.size), dtype=float)
    outdat = utils.add_col(outdat, 'z0', np.array([0]*outdat.size), dtype=float)

    # Now start the outer loop, over halos. Sorry.
    i0=0
    for i, g in enumerate(groups):

    	# MPI handling
        if i%size!=rank:
            continue

        # Select the galaxies in this particular halo
        select = (group_ids==g)
        ngal = len(data[select])
        if verbose:
        	print 'Halo %g contains %d object(s)'%(g,ngal)

        outdat['halo_id'][select] = group_ids[select]

        # We might want to apply some threshold, and ignore halos with fewer galaxies than some threshold 
        if (ngal<nmin):
            if verbose:
            	print 'Skipping halo %d as it contains <%d galaxies'%(g,nmin)
            continue

        # Now apply the symmetrisation operation to these galaxies
        symmetrised_halo = symmetrise_halo5(data[select], verbose=True, g=g)

        if verbose:
            print g, ngal

        # Transfer the new symmetrised columns to the output array
        for name in outdat.dtype.names:
            outdat[name][select] = symmetrised_halo[name]

        # Count the number of halos actually altered by the process
        i0+=ngal

    # Finally store the result as a FITS file
    outfits = fi.FITS(savedir.replace('.fits', '%d.fits'%rank), 'rw')
    outfits.write(outdat)
    outfits.close()

    print 'Done'
    return 0

def symmetrise_halo4(data, verbose=True, g=None):
    """
	Takes an array of objects associated with a particular halo,
	queries the database on coma to find the halo centroid positions
	and then applies a random rotation to each.
	The output is an array in the same format as given, but with
	new positions and shapes."""

    N = data['npart_baryon']
    n = data.size

    # Work out the centroid about which to rotate  
    data, positions, (x0, y0, z0) = get_wrapped_positions(g, data) 
    cent = np.array([x0,y0,z0])

    # New array to store the output
    rot = np.zeros(data.size, dtype=data.dtype)  
    for name in data.dtype.names:
        rot[name][:n] = data[name]

    rot = utils.add_col(rot, 'x0', np.array([x0]*data.size), dtype=float)
    rot = utils.add_col(rot, 'y0', np.array([y0]*data.size), dtype=float)
    rot = utils.add_col(rot, 'z0', np.array([z0]*data.size), dtype=float)

    # Start the inner loop, over individual galaxies. Sorry.
    for i in xrange(n):

    	# Cartesian galaxy position, relative to the halo centre
    	pos = positions.T[i]

    	# Galaxy position in spherical polar coordinates
    	R = np.sqrt(sum(pos*pos))
    	phi = np.arccos(pos[2]/R) * pos[0]/abs(pos[0])
    	theta = np.arcsin(pos[1]/R/np.sin(phi))
    	rot['rh'][i] = R

        # We can skip the rest for central objects
        if (data['most_massive'][i]==1):
            rot = check_wrapping(i, rot)
            if verbose:
                print 'skipping object -- it is classified as a central galaxy'
            continue

        if verbose:
            print i,

        # Construct a 3x3 matrix to rotate galaxies by a random angle in the range [0, 2/pi]
        # about a randomly chosen axis.
        Rxyz = build_rotation_matrix()

        # Rotated positions
        rotated = np.dot(Rxyz,pos)
        rot['x'][i] = copy.deepcopy(rotated[0])+x0
        rot['y'][i] = copy.deepcopy(rotated[1])+y0
        rot['z'][i] = copy.deepcopy(rotated[2])+z0

        if verbose:
            print 'New position (x, y, z) : %3.3f, %3.3f %3.3f'%(rot['x'][i],rot['y'][i],rot['z'][i])

        # Apply the same rotation to the three orientation vectors
        a3d = np.array([data['a1'][i], data['a2'][i], data['a3'][i]])
        b3d = np.array([data['b1'][i], data['b2'][i], data['b3'][i]])
        c3d = np.array([data['c1'][i], data['c2'][i], data['c3'][i]])
        arot = np.dot(Rxyz,a3d)
        brot = np.dot(Rxyz,b3d)
        crot = np.dot(Rxyz,c3d)
        axes = {'a':arot,'b':brot,'c':crot}

        # And for the dark matter component
        a3d_dm = np.array([data['a1_dm'][i], data['a2_dm'][i], data['a3_dm'][i]])
        b3d_dm = np.array([data['b1_dm'][i], data['b2_dm'][i], data['b3_dm'][i]])
        c3d_dm = np.array([data['c1_dm'][i], data['c2_dm'][i], data['c3_dm'][i]])
        arot_dm = np.dot(Rxyz,a3d_dm)
        brot_dm = np.dot(Rxyz,b3d_dm)
        crot_dm = np.dot(Rxyz,c3d_dm)
        axes_dm = {'a':arot_dm,'b':brot_dm,'c':crot_dm}

        # Stays the same
        rot['npart_baryon'][i] = data['npart_baryon'][i]

        # We also need to shift a few objects near the edges, which the rotation leaves outside the simulation box,
        # back to the other side of the universe 
        if verbose: 
            print 'Rotation leaves %d object(s) outside the simulation box.'%(rot['x'][(rot['x']<0) | (rot['y']<0) | (rot['z']<0)].size)

        rot = check_wrapping(i, rot)

        if rot['x'][i]<0 : import pdb ; pdb.set_trace()

        # Store everything in the appropriate place
        for j in xrange(3):
        	for axis in ['a','b','c']:
        		rot['%c%d'%(axis,j+1)][i] =  axes[axis][j]
        		rot['%c%d_dm'%(axis,j+1)][i] =  axes_dm[axis][j]

        # Work out the projected 2D ellipticities for the stellar component
        rot = project_ellipticities(i, rot, suffix='')
        # Same for dark matter 
        rot = project_ellipticities(i, rot, suffix='dm')
        
    return rot


def symmetrise_halo5(data, verbose=True, g=None):
    """
    Takes an array of objects associated with a particular halo,
    queries the database on coma to find the halo centroid positions
    and then applies a random rotation to each.
    The output is an array in the same format as given, but with
    new positions and shapes."""

    N = data['npart_baryon']
    n = data.size

    # Work out the centroid about which to rotate  
    data, positions, (x0, y0, z0) = get_wrapped_positions(g, data) 
    cent = np.array([x0,y0,z0])

    # New array to store the output
    rot = np.zeros(data.size, dtype=data.dtype)  
    for name in data.dtype.names:
        rot[name][:n] = data[name]

    rot = utils.add_col(rot, 'x0', np.array([x0]*data.size), dtype=float)
    rot = utils.add_col(rot, 'y0', np.array([y0]*data.size), dtype=float)
    rot = utils.add_col(rot, 'z0', np.array([z0]*data.size), dtype=float)

    # Start the inner loop, over individual galaxies. Sorry.
    for i in xrange(n):

        # Cartesian galaxy position, relative to the halo centre
        pos = positions.T[i]

        # Galaxy position in spherical polar coordinates
        R = np.sqrt(sum(pos*pos))
        phi = np.arccos(pos[2]/R) * pos[0]/abs(pos[0])
        theta = np.arcsin(pos[1]/R/np.sin(phi))
        rot['rh'][i] = R

        # We can skip the rest for central objects
        if (data['most_massive'][i]==1):
            rot = check_wrapping(i, rot)
            if verbose:
                print 'skipping object -- it is classified as a central galaxy'
            continue

        if verbose:
            print i,

        # Choose a random point on a sphere about the centroid of radius R
        rotated = sample_sphere(1, norm=R)
        rotated = np.array([rotated[0][0], rotated[1][0], rotated[2][0]])
        rot['x'] = rotated[0]
        rot['y'] = rotated[1]
        rot['z'] = rotated[2]

        # Then work backwards to get the rotation matrix needed to transform the initial position to the rotated one
        rotation_axis, rotation_angle = infer_rotation_angle(pos,rotated)
        Rxyz = build_rotation_matrix(alpha=rotation_angle, vec=rotation_axis)

        if verbose:
            print 'New position (x, y, z) : %3.3f, %3.3f %3.3f'%(rot['x'][i],rot['y'][i],rot['z'][i])

        # Apply the same rotation to the three orientation vectors
        a3d = np.array([data['a1'][i], data['a2'][i], data['a3'][i]])
        b3d = np.array([data['b1'][i], data['b2'][i], data['b3'][i]])
        c3d = np.array([data['c1'][i], data['c2'][i], data['c3'][i]])
        arot = np.dot(Rxyz,a3d)
        brot = np.dot(Rxyz,b3d)
        crot = np.dot(Rxyz,c3d)
        axes = {'a':arot,'b':brot,'c':crot}

        # And for the dark matter component
        a3d_dm = np.array([data['a1_dm'][i], data['a2_dm'][i], data['a3_dm'][i]])
        b3d_dm = np.array([data['b1_dm'][i], data['b2_dm'][i], data['b3_dm'][i]])
        c3d_dm = np.array([data['c1_dm'][i], data['c2_dm'][i], data['c3_dm'][i]])
        arot_dm = np.dot(Rxyz,a3d_dm)
        brot_dm = np.dot(Rxyz,b3d_dm)
        crot_dm = np.dot(Rxyz,c3d_dm)
        axes_dm = {'a':arot_dm,'b':brot_dm,'c':crot_dm}

        # Stays the same
        rot['npart_baryon'][i] = data['npart_baryon'][i]

        # We also need to shift a few objects near the edges, which the rotation leaves outside the simulation box,
        # back to the other side of the universe 
        if verbose: 
            print 'Rotation leaves %d object(s) outside the simulation box.'%(rot['x'][(rot['x']<0) | (rot['y']<0) | (rot['z']<0)].size)

        rot = check_wrapping(i, rot)

        if rot['x'][i]<0 : import pdb ; pdb.set_trace()

        # Store everything in the appropriate place
        for j in xrange(3):
            for axis in ['a','b','c']:
                rot['%c%d'%(axis,j+1)][i] =  axes[axis][j]
                rot['%c%d_dm'%(axis,j+1)][i] =  axes_dm[axis][j]

        # Work out the projected 2D ellipticities for the stellar component
        rot = project_ellipticities(i, rot, suffix='')
        # Same for dark matter 
        rot = project_ellipticities(i, rot, suffix='dm')
        
    return rot

def project_ellipticities(i, rot, suffix=''):

	if (suffix!=''):
		suffix = '_%s'%suffix

	a3d = np.array([[rot['a1'+suffix][i], rot['a2'+suffix][i], rot['a3'+suffix][i]]])
	b3d = np.array([[rot['b1'+suffix][i], rot['b2'+suffix][i], rot['b3'+suffix][i]]])
	c3d = np.array([[rot['c1'+suffix][i], rot['c2'+suffix][i], rot['c3'+suffix][i]]])
	q3d = np.array([np.sqrt(rot['lambda_a'+suffix][i]/rot['lambda_c'+suffix][i])])
	s3d = np.array([np.sqrt(rot['lambda_b'+suffix][i]/rot['lambda_c'+suffix][i])])

	e1,e2 = utils.project_3d_shape(a3d, b3d, c3d, q3d, s3d)
	rot['e1'+suffix][i] = e1
	rot['e2'+suffix][i] = e2

	return rot


def get_wrapped_positions(g, data, use_database=True):
    # Query the DB for the halo centre
    if use_database:
    	info = find_centre(g)
    else:
    	info = np.zeros(1,dtype=['x','y','z'])
    	mask = (data['most_massive']==1)
    	info['x'] = data['x'][mask]
    	info['y'] = data['y'][mask]
    	info['z'] = data['z'][mask]
    
    for comp in ['x','y','z']:
        if (data[comp]<5).max() and (data[comp]>95).max():
            # Shift everything back down to the lower edge of the simulation box
            # This will produce negative positions, which is fine if we shift the
            # galaxies outside the box back after the rotation

            # The 5 Mpc border is a bit arbitary, but it's probably safe to say any
            # halo with galaxies within both the uppermost and lowermost 5 Mpc is
            # an artefact of the periodic boundary conditions    
            data[comp][data[comp]>95]-=100

        # This should handle the case where the centroid is on the opposite
        # side of the box to all of the galaxies surviving cuts 
        if ( ((info[comp]/1000)>95) and (data[comp]<5).max() ) or ( ((info[comp]/1000)<5) and (data[comp]>95).max() ):
                info[comp]-=100*1000

    x0,y0,z0 = info['x']/1000., info['y']/1000., info['z']/1000.
    positions = np.array([data['x']-x0, data['y']-y0, data['z']-z0])

    return data, positions, (x0,y0,z0)

def sample_sphere(npoints, norm=1, seed=9000):
    """Generates random points on a sphere of radius R"""
    np.random.seed(seed)
    u1 = np.random.rand(npoints)
    u2 = np.random.rand(npoints)

    # Define a vector position
    theta = np.arccos(2*u1-1)
    phi = 2 * np.pi * u2

    x = np.sin(theta)*np.cos(phi)
    y = np.sin(theta)*np.sin(phi)
    z = np.cos(theta)
    vec= np.array([x,y,z])
    vec *= norm

    return vec

def infer_rotation_angle(position0, position):
    """Work out two rotation angles that will transform one given 
    3D vector into another."""

    # Work out the unit vector orthogonal to the plane defined by the two position vectors
    rotation_axis = np.cross(position0, position)
    rotation_axis /= utils.get_norm(rotation_axis)

    # Then calculate the angle between the two position vectors in the 2D plane defined by them
    alpha = np.arccos( np.dot(position/utils.get_norm(position), position0/utils.get_norm(position0)) )
    cross = np.cross(position, position0);

    #alpha = -alpha

    return rotation_axis, alpha



def build_rotation_matrix(theta=None, phi=None, alpha=None, vec=[]):
    """Generates a random rotation matrix, 
    which transforms a given 3D position vector to a random position on the sphere.
    Usage : rotated = R.unrotated"""
    u1 = np.random.rand()
    u2 = np.random.rand()

    # Use two random values from a uniform distribution to generate an axis about which to rotate
    if phi==None:
        phi = np.arccos(2*u1-1)
    if theta==None:
        theta = 2 * np.pi * u2

    # Work out a Cartesian unit vector defined by the rotation axis 
    if (len(vec)==0):
        x = np.sin(theta)*np.cos(phi)
        y = np.sin(theta)*np.sin(phi)
        z = np.cos(theta)
        vec= np.array([x,y,z])
    vec/=utils.get_norm(vec)

    # Now generate a random rotation angle about that axis
    if alpha==None:
        alpha = np.random.rand() * 2 * np.pi
    cosa = np.cos(alpha)
    sina = np.sin(alpha)
    ux = vec[0]
    uy = vec[1]
    uz = vec[2]

    # Construct the 3x3 rotation matrix
    R = np.array([ [ (cosa + ux*ux*(1-cosa)), (ux*uy*(1-cosa) - uz*sina), (ux*uz*(1-cosa) + uy*sina)],
    [uy*ux*(1-cosa) + uz*sina, cosa + uy*uy*(1-cosa), uy*uz*(1-cosa) - ux*sina],
    [uz*ux*(1-cosa) - uy*sina, uz*uy*(1-cosa) + ux*sina, cosa + uz*uz*(1-cosa) ]] )

    return R


def find_centre(g):
    """Locate the mass centroid of a particular group of galaxies"""
    import pymysql as mdb
    sql = 'SELECT x,y,z FROM subfind_groups WHERE groupId=%d AND snapnum=85;'%g

    sqlserver='localhost' ; user='flanusse' ; password='mysqlpass' ; dbname='mb2_hydro'
    unix_socket='/home/rmandelb.proj/flanusse/mysql/mysql.sock'
    db = mdb.connect(sqlserver, user, password, dbname, unix_socket=unix_socket)

    c = db.cursor()
    c.execute(sql)
    results = fromarrays(np.array(c.fetchall()).squeeze().T,names='x,y,z')
    return results

def check_wrapping(i, rot, boxsize=100):
	# Lower edge on each axis
	if (rot['x'][i]<0):
		rot['x'][i]+=boxsize
	if (rot['y'][i]<0):
		rot['y'][i]+=boxsize
	if (rot['z'][i]<0):
		rot['z'][i]+=boxsize

	# Upper edge on each axis
	if (rot['x'][i]>boxsize):
		rot['x'][i]-=boxsize
	if (rot['y'][i]>boxsize):
		rot['y'][i]-=boxsize
	if (rot['z'][i]>boxsize):
		rot['z'][i]-=boxsize

	return rot


