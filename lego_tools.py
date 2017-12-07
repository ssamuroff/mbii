import halotools as ht
import numpy as np
import copy
import pylab as plt
plt.switch_backend('agg')
import halotools.mock_observables as pretending
from mbii.readsubhalo import *
from mbii.properties import *


import numpy as np
#import astropy.table as tb
#import astropy.io.fits as pf


def plot_ellipse(semimaj=1,semimin=1,phi=0,x_cent=0,y_cent=0,theta_num=1e3,ax=None,plot_kwargs=None,\
                    fill=False,fill_kwargs=None,data_out=False,cov=None,mass_level=0.68):
    '''
        The function creates a 2D ellipse in polar coordinates then transforms to cartesian coordinates.
        It can take a covariance matrix and plot contours from it.
        
        semimaj : float
            length of semimajor axis (always taken to be some phi (-90<phi<90 deg) from positive x-axis!)

        semimin : float
            length of semiminor axis

        phi : float
            angle in radians of semimajor axis above positive x axis

        x_cent : float
            X coordinate center

        y_cent : float
            Y coordinate center

        theta_num : int
            Number of points to sample along ellipse from 0-2pi

        ax : matplotlib axis property
            A pre-created matplotlib axis

        plot_kwargs : dictionary
            matplotlib.plot() keyword arguments

        fill : bool
            A flag to fill the inside of the ellipse 

        fill_kwargs : dictionary
            Keyword arguments for matplotlib.fill()

        data_out : bool
            A flag to return the ellipse samples without plotting

        cov : ndarray of shape (2,2)
            A 2x2 covariance matrix, if given this will overwrite semimaj, semimin and phi

        mass_level : float
            if supplied cov, mass_level is the contour defining fractional probability mass enclosed
            for example: mass_level = 0.68 is the standard 68% mass

    '''
    # Get Ellipse Properties from cov matrix
    if cov is not None:
        eig_vec,eig_val,u = np.linalg.svd(cov)
        # Make sure 0th eigenvector has positive x-coordinate
        if eig_vec[0][0] < 0:
            eig_vec[0] *= -1
        semimaj = np.sqrt(eig_val[0])
        semimin = np.sqrt(eig_val[1])
        if mass_level is None:
            multiplier = np.sqrt(2.279)
        else:
            distances = np.linspace(0,20,20001)
            chi2_cdf = chi2.cdf(distances,df=2)
            multiplier = np.sqrt(distances[np.where(np.abs(chi2_cdf-mass_level)==np.abs(chi2_cdf-mass_level).min())[0][0]])
        semimaj *= multiplier
        semimin *= multiplier
        phi = np.arccos(np.dot(eig_vec[0],np.array([1,0])))
        if eig_vec[0][1] < 0 and phi > 0:
            phi *= -1

    # Generate data for ellipse structure
    theta = np.linspace(0,2*np.pi,theta_num)
    r = 1 / np.sqrt((np.cos(theta))**2 + (np.sin(theta))**2)
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    data = np.array([x,y])
    S = np.array([[semimaj,0],[0,semimin]])
    R = np.array([[np.cos(phi),-np.sin(phi)],[np.sin(phi),np.cos(phi)]])
    T = np.dot(R,S)
    data = np.dot(T,data)
    data[0] += x_cent
    data[1] += y_cent

    # Output data?
    if data_out == True:
        return data

    # Plot!
    return_fig = False
    if ax is None:
        return_fig = True
        fig,ax = plt.subplots()

    if plot_kwargs is None:
        ax.plot(data[0],data[1],color='purple',linestyle='-')
    else:
        ax.plot(data[0],data[1],**plot_kwargs)

    if fill == True:
        ax.fill(data[0],data[1],**fill_kwargs)

    if return_fig == True:
        return fig


def symmetrise_halo(data, nrot=100, gids=None):
    print 'Will apply %d random rotations'%nrot
    np.random.seed(9000)

    # Work out the centroid about which to rotate
    
    x0 = np.mean(data['x'][np.isfinite(data['x'])])
    y0 = np.mean(data['y'][np.isfinite(data['y'])])
    z0 = np.mean(data['z'][np.isfinite(data['z'])])
    cent = np.array([x0,y0,z0])
    positions = np.array([data['x']-x0, data['y']-y0, data['z']-z0])
    c = np.array([data['c3'], data['c2'], data['c1']])
    a = np.array([data['a1'], data['a2'], data['a3']])
    b = np.array([data['b3'], data['b2'], data['b1']])

    rot = np.zeros(data.size*(nrot+1), dtype=data.dtype)
    n = data.size
    for name in data.dtype.names:
        rot[name][:n] = data[name]


    for i in xrange(nrot):
        print i,
        # Define a random rotation about each axis
        alpha = 2 * np.pi * (np.random.rand() - 0.5) 
        beta = 2 * np.pi * (np.random.rand() - 0.5) 
        gamma = 2 * np.pi * (np.random.rand() - 0.5)

        print alpha, beta, gamma
        #alpha,beta,gamma=np.pi,0,0

        #import pdb ; pdb.set_trace()

        # Rotate the position vector
        Rx = np.array([[1,0,0], [0, np.cos(alpha), -np.sin(alpha)], [0, np.sin(alpha), np.cos(alpha)]])
        rotated = np.dot(Rx,positions)
        Ry = np.array([[np.cos(beta), 0, np.sin(beta)],[0, 1, 0],[-np.sin(beta), 0, np.cos(beta)]])
        rotated = np.dot(Ry,rotated)
        Rz = np.array([[np.cos(gamma), -np.sin(gamma), 0], [np.sin(gamma), np.cos(gamma), 0], [0,0,1]])
        rotated = np.dot(Rz,rotated)

        # Now do the shapes
        arot = np.dot(Rx,a)
        arot = np.dot(Ry,arot)
        arot = np.dot(Rz,arot)

        brot = np.dot(Rx,b)
        brot = np.dot(Ry,brot)
        brot = np.dot(Rz,brot)

        crot = np.dot(Rx,c)
        crot = np.dot(Ry,crot)
        crot = np.dot(Rz,crot)

        rot['npart'][(i+1)*n:(i+2)*n] = data['npart']

        
        rot['x'][(i+1)*n:(i+2)*n] = copy.deepcopy(rotated[0])+x0
        rot['y'][(i+1)*n:(i+2)*n] = copy.deepcopy(rotated[1])+y0
        rot['z'][(i+1)*n:(i+2)*n] = copy.deepcopy(rotated[2])+z0

        rot['c3'][(i+1)*n:(i+2)*n]=crot[0]
        rot['c2'][(i+1)*n:(i+2)*n]=crot[1]
        rot['c1'][(i+1)*n:(i+2)*n]=crot[2]
        rot['b3'][(i+1)*n:(i+2)*n]=brot[0]
        rot['b2'][(i+1)*n:(i+2)*n]=brot[1]
        rot['b1'][(i+1)*n:(i+2)*n]=brot[2]
        rot['a3'][(i+1)*n:(i+2)*n]=arot[0]
        rot['a2'][(i+1)*n:(i+2)*n]=arot[1]
        rot['a1'][(i+1)*n:(i+2)*n]=arot[2]

    import pdb ; pdb.set_trace()

    return rot


def visualise_halo(data, group=0, mask=None, h=None, axes=['x','y'], cmap='Purples', shapes=False, projected=None, gid=None, axis=True, ticks=False, fontsize=8):

    # Load the subhalo data if missing
    if h is None:
        root_folder='/physics/yfeng1/mb2'
        snapshot='085'
        snap = SnapDir(snapshot, root_folder)
        h = snap.readsubhalo()

    if mask is None:
        print 'No mask provided. Using all objects'
        mask = np.ones(data.size).astype(bool)

    if gid is None:
        gid = np.array(h['groupid'])

    if axis: 
        axis = plt.subplot(111,aspect='equal')

    if not shapes:
        plt.scatter(data[axes[0]][mask & (gid==group)], data[axes[1]][mask & (gid==group)], marker='.', cmap='Purples', c=np.log10(data['npart'][mask & (gid==group)]*2.2e6)  )
    else:
        angle = np.arctan2(projected['a1'],projected['a2'])[mask & (gid==group)]

        for (semimaj, semimin, phi, x0, y0) in zip(projected['lambda2'][mask & (gid==group)],projected['lambda1'][mask & (gid==group)], angle, data[axes[0]][mask & (gid==group)], data[axes[1]][mask & (gid==group)]):
            #import pdb ; pdb.set_trace()
            q = semimin/semimaj
            plot_ellipse(semimaj=0.005e4/q,semimin=0.005e4*q,phi=phi,x_cent=x0,y_cent=y0, ax=axis)


    plt.yticks(fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    if not ticks:
        plt.yticks(visible=False)
        plt.xticks(visible=False)
    else:
        plt.xlabel('$%s$ / $h^{-1}$ kpc'%axes[0], fontsize=fontsize)
        plt.ylabel('$%s$ / $h^{-1}$ kpc'%axes[1], fontsize=fontsize) 


    dx = 8500.0
    plt.xlim(data[axes[0]][mask & (gid==group)].mean() - dx/2, data[axes[0]][mask & (gid==group)].mean() + dx/2)
    plt.ylim(data[axes[1]][mask & (gid==group)].mean() - dx/2, data[axes[1]][mask & (gid==group)].mean() + dx/2)

    if axis:
        plt.subplots_adjust(top=0.98,bottom=0.14, right=0.98, left=0.14)



def cartesian_to_cylindrical(vec):
    theta = np.arctan2(vec[1],vec[0])
    r = vec[1]*np.sin(theta) - vec[0]*np.cos(theta)
    t = vec[1]*np.cos(theta) - vec[0]*np.sin(theta)
    return r, t, vec[2]


def process_bootstrap(bootstrap_iterations):
    mu = np.array(bootstrap_iterations).mean(axis=0)
    res = (np.array(bootstrap_iterations) - mu)**2
    niter = len(bootstrap_iterations)

    B = np.sqrt(sum(res)/(niter-1))

    return B

def add_col(rec, name, arr=[], dtype=None, verbose=False):
    """Generic function to add a new column to a structured numpy array."""

    if name in rec.dtype.names:
        if verbose:
            print "Table already has a column called %s"%name
        return rec

    if len(arr)==0:
        arr=np.zeros(len(rec))

    arr = np.asarray(arr)
    if dtype is None:
        dtype = arr.dtype

    newdtype = np.dtype(rec.dtype.descr + [(name, dtype)])
    newrec = np.empty(rec.shape, dtype=newdtype)
    for field in rec.dtype.fields:
        newrec[field] = rec[field]

    newrec[name] = arr

    return newrec



def construct_random_cat(data, mask=None, format='halotools', f=1):
	"""Construct a catalogue of random points drawn from the same volume as the data provided."""

	if f!=1:
		print 'Random point multiplier: %3.3f'%f
	if mask is None:
		mask = np.ones(data.size).astype(bool)
	rx = (np.random.random(size=int(data['x'][mask].size*f)) - 0.5) * (data['x'][mask].max()-data['x'][mask].min()) + data['x'][mask].mean()
	ry = (np.random.random(size=int(data['x'][mask].size*f)) - 0.5) * (data['y'][mask].max()-data['y'][mask].min()) + data['y'][mask].mean()
	rz = (np.random.random(size=int(data['x'][mask].size*f)) - 0.5) * (data['z'][mask].max()-data['z'][mask].min()) + data['z'][mask].mean()

	if format.lower()=='halotools':
		return pretending.return_xyz_formatted_array(rx, ry, rz)

	elif format.lower()=='treecorr':
		return treecorr.Catalog(x=rx, y=ry, z=rz)


def choose_cs_mask(data, ctype):

    if ctype=='c':
        return (data['central']==1)
    if ctype=='s':
        return (data['central']!=1)
    if ctype=='a':
        return np.ones(data['central'].size).astype(bool)

    print 'Unrecognised galaxy type:', ctype
    return None



def nfw(r,c,a):
    return a/(r/c)/(1+r/c)/(1+r/c)


def convert_to_halotools(data, snapshot=85, read_mass=False):
	z = snapshots[snapshot]
	names = data.dtype.names
	if ('mass' in names) and read_mass:
		mass = data['mass']
	else:
		mass =0

	catalogue = ht.sim_manager.UserSuppliedPtclCatalog(redshift=z, Lbox=100, particle_mass=mass, x=data['x']/1e3, y=data['y']/1e3, z=data['z']/1e3, ptcl_ids=data['subfindId'])

	return catalogue


def parse_mask(mask,array):
    if mask is None:
        return np.ones(array.size).astype(bool)
    else:
        return mask

def footprint_sub(ra,dec,rasep,decsep,nside,fig, cmap='none'):
    plt.close() ; fig = plt.figure(figsize=(6.5,6))

    import skymapper as skm

    bc, ra0, dec0, vertices = skm.getCountAtLocations(ra, dec, nside=nside, return_vertices=True)

    # setup figure
    import matplotlib.cm as cm
    if cmap=='none':
        cmap = cm.YlOrRd
    ax = fig.add_subplot(111, aspect='equal')

    # setup map: define AEA map optimal for given RA/Dec
    proj = skm.createConicMap(ax, ra0, dec0, proj_class=skm.AlbersEqualAreaProjection)
    # add lines and labels for meridians/parallels (separation 5 deg)
    sep = 5
    meridians = np.arange(-90, 90+decsep, decsep)
    parallels = np.arange(0, 360+rasep, rasep)
    skm.setMeridianPatches(ax, proj, meridians, linestyle='-', lw=0.5, alpha=0.3, zorder=2)
    skm.setParallelPatches(ax, proj, parallels, linestyle='-', lw=0.5, alpha=0.3, zorder=2)
    skm.setMeridianLabels(ax, proj, meridians, loc="left", fmt=skm.pmDegFormatter)
    skm.setParallelLabels(ax, proj, parallels, loc="top")

    # add vertices as polygons
    vmin, vmax = np.percentile(bc,[10,90])
    poly = skm.addPolygons(vertices, proj, ax, color=bc, vmin=vmin, vmax=vmax, cmap=cmap, zorder=3, rasterized=True)

    # add colorbar
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.0)
    cb = fig.colorbar(poly, cax=cax)
    cb.set_label('$n_g$ [arcmin$^{-2}$]')
    cb.solids.set_edgecolor("face")

    skm.addFootprint('DES', proj, ax, zorder=10, edgecolor='#2222B2', facecolor='None', lw=2)

    #Fiddle with axis limit so we can see the whole range
    xmin,_=proj(100,-30)
    _,xmax=ax.get_xlim()
    ax.set_xlim(xmin,xmax)
    ymin,ymax=ax.get_ylim()
    r=ymax-ymin
    ymin-=r/10.
    ymax+=r/5.
    ax.set_ylim(ymin,ymax)
    return

def project_to_sky(data):
    r = np.sqrt(data['x']*data['x'] + data['y']*data['y'] + data['z']*data['z'])
    ra = np.arccos(data['z']/r)
    dec = np.arccos(data['x']/r/np.sin(ra))

    ra = ra * 180./np.pi
    dec = dec * 180./np.pi

    fig = plt.figure(0)
    footprint_sub(ra, dec, 10, 10, 2048, fig, cmap='Purples')

    return fig

