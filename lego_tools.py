import halotools as ht
import numpy as np
import halotools.mock_observables as pretending

import numpy as np
#import astropy.table as tb
#import astropy.io.fits as pf

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

