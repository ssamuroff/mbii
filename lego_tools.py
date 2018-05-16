import halotools as ht
import numpy as np
import copy
import fitsio as fi
import pylab as plt
plt.switch_backend('agg')
#import halotools.mock_observables as pretending
from mbii.readsubhalo import *
from mbii.properties import *
#import mbii.tests as tests
import glob
import numpy as np
import treecorr
#import astropy.table as tb
#import astropy.io.fits as pf
from numpy.core.records import fromarrays
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def visualise_3d(x,y,z):  
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, '.', color='purple')


def equalise_binning(data1, data2, rmin, rmax, nbin, tol=1000.):
    """Iteratively shift the bin edges and redo the pair counting until the bins
    contain equal numbers of pairs."""

    from halotools.mock_observables.pair_counters import npairs_3d

    rbins = np.logspace(np.log10(rmin), np.log10(rmax), nbin)
    upper = rbins[1:]
    lower = rbins[:-1]

    positions1 = np.vstack((data1['x'],data1['y'],data1['z'])).T
    positions2 = np.vstack((data2['x'],data2['y'],data2['z'])).T
    counts = npairs_3d(positions1, positions2, rbins,
                       verbose=False, num_threads=1,
                       approx_cell1_size=None,
                       approx_cell2_size=None)

    nmax = counts.max()
    nmin = counts.min()
    R = nmax*1./nmin
    iiter=0

    while R>tol:
        rnew=[]
        imax=np.argwhere(counts==counts.max())[0,0]
        imin=np.argwhere(counts==counts.min())[0,0]
        for j,r in enumerate(rbins):
            if j-1==imax-1:
                r0 = np.sqrt(rbins[j-1]*rbins[j])
                rnew += [r0,r]
            elif j-1==imin:
                continue
            else:
                rnew.append(r)

        counts = npairs_3d(positions1, positions2, rnew,
                    verbose=False, num_threads=1,
                    approx_cell1_size=None,
                    approx_cell2_size=None)

        nmax = counts.max()
        nmin = counts.min()
        R = nmax*1./nmin
        iiter+=1
        rbins = rnew
        print 'Done iteration ', iiter, R

    return np.array(rbins)
    

def find_equal_bin_edges(x,nbins,w=None):
    """
    For an array x, returns the boundaries of nbins equal (possibly weighted by w) bins.
    """

    if w is None:
      xs=np.sort(x)
      r=np.linspace(0.,1.,nbins+1.)*(len(x)-1)
      return xs[r.astype(int)]

    fail=False
    ww=np.sum(w)/nbins
    i=np.argsort(x)
    k=np.linspace(0.,1.,nbins+1.)*(len(x)-1)
    k=k.astype(int)
    r=np.zeros((nbins+1))
    ist=0
    for j in xrange(1,nbins):
      # print k[j],r[j-1]
      if k[j]<r[j-1]:
        print 'Random weight approx. failed - attempting brute force approach'
        fail=True
        break
      w0=np.sum(w[i[ist:k[j]]])
      if w0<=ww:
        for l in xrange(k[j],len(x)):
          w0+=w[i[l]]
          if w0>ww:
            r[j]=x[i[l]]
            ist=l
            break
      else:
        for l in xrange(k[j],0,-1):
          w0-=w[i[l]]
          if w0<ww:
            r[j]=x[i[l]]
            ist=l
            break

    if fail:

      ist=np.zeros((nbins+1))
      ist[0]=0
      for j in xrange(1,nbins):
        wsum=0.
        for k in xrange(ist[j-1].astype(int),len(x)):
          wsum+=w[i[k]]
          if wsum>ww:
            r[j]=x[i[k-1]]
            ist[j]=k
            break

    r[0]=x[i[0]]
    r[-1]=x[i[-1]]

    return equal_r

    

def find_equal_bin_edges(x,nbins,w=None):
    """
    For an array x, returns the boundaries of nbins equal (possibly weighted by w) bins.
    """

    if w is None:
      xs=np.sort(x)
      r=np.linspace(0.,1.,nbins+1.)*(len(x)-1)
      return xs[r.astype(int)]

    fail=False
    ww=np.sum(w)/nbins
    i=np.argsort(x)
    k=np.linspace(0.,1.,nbins+1.)*(len(x)-1)
    k=k.astype(int)
    r=np.zeros((nbins+1))
    ist=0
    for j in xrange(1,nbins):
      # print k[j],r[j-1]
      if k[j]<r[j-1]:
        print 'Random weight approx. failed - attempting brute force approach'
        fail=True
        break
      w0=np.sum(w[i[ist:k[j]]])
      if w0<=ww:
        for l in xrange(k[j],len(x)):
          w0+=w[i[l]]
          if w0>ww:
            r[j]=x[i[l]]
            ist=l
            break
      else:
        for l in xrange(k[j],0,-1):
          w0-=w[i[l]]
          if w0<ww:
            r[j]=x[i[l]]
            ist=l
            break

    if fail:

      ist=np.zeros((nbins+1))
      ist[0]=0
      for j in xrange(1,nbins):
        wsum=0.
        for k in xrange(ist[j-1].astype(int),len(x)):
          wsum+=w[i[k]]
          if wsum>ww:
            r[j]=x[i[k-1]]
            ist[j]=k
            break

    r[0]=x[i[0]]
    r[-1]=x[i[-1]]

    return r



def symmetrise_catalogue3(data=None, seed=4000, mask=None,filename='/home/ssamurof/massive_black_ii/subhalo_cat-nthreshold5.fits', savedir=None, replace=False, verbose=True, nmin=0, rank=0, size=1):
    root_folder='/physics/yfeng1/mb2'
    snapshot='085'

    np.random.seed(99000)

    if mask is None:
        mask = np.ones(data.size).astype(bool)

    print savedir
    if savedir==filename:
        print "WARNING: target and input file paths are the same"
        import pdb ; pdb.set_trace()

    group_ids = data['halo_id'][mask]
    groups = np.unique(group_ids)


    # Copy over the old data
    outdat = np.zeros(len(data), dtype=data.dtype)
    for name in outdat.dtype.names:
        outdat[name] = data[name]

    outdat = add_col(outdat, 'x0', np.array([0]*outdat.size), dtype=float)
    outdat = add_col(outdat, 'y0', np.array([0]*outdat.size), dtype=float)
    outdat = add_col(outdat, 'z0', np.array([0]*outdat.size), dtype=float)

    i0=0
    for i, g in enumerate(groups):

        if i%size!=rank:
            continue
        
        select = (group_ids==g)

        ngal = len(data[select])
        outdat['halo_id'][select] = group_ids[select]
        if ngal<nmin:
            print 'Skipping halo %d as it contains <%d galaxies'%(g,nmin)
            continue

        symmetrised_halo = symmetrise_halo4(data[select], gids=group_ids[select], verbose=True, g=g)

        if verbose:
            print g, ngal

        for name in outdat.dtype.names:
            outdat[name][select] = symmetrised_halo[name]

        i0+=ngal

    outfits = fi.FITS(savedir.replace('.fits', '%d.fits'%rank), 'rw')
    outfits.write(outdat)
    outfits.close()

    print 'Done'
    return 0

def get_wrapped_positions(g, data):
    # Query the DB for the halo centre
    info = find_centre(g)
    
    for comp in ['x','y','z']:
        if (data[comp]<5).max() and (data[comp]>95).max():
            # Shift everything back down to the lower edge of the simulation box
            # This will produce negative positions, which is fine if we shift the
            # galaxies outside the box back after the rotation
            data[comp][data[comp]>95]-=100

        # This should handle the case where the centroid is on the opposite
        # side of the box to all of the galaxies surviving cuts 
        if ( ((info[comp]/1000)>95) and (data[comp]<5).max() ) or ( ((info[comp]/1000)<5) and (data[comp]>95).max() ):
                info[comp]-=100*1000

    x0,y0,z0 = info['x']/1000., info['y']/1000., info['z']/1000.
    positions = np.array([data['x']-x0, data['y']-y0, data['z']-z0])

    return data, positions, (x0,y0,z0)

def symmetrise_halo4(data, gids=None, verbose=True, g=None):
    if verbose:
        print 'Will apply one random rotation per subhalo'

    # Work out the centroid about which to rotate
    import weightedstats as ws
    N = data['npart_baryon']
    # Let's try something slightly different here
    # get the centroid position from the database rather than trying to recalculate it
    
    data, positions, (x0, y0, z0) = get_wrapped_positions(g, data) 
    cent = np.array([x0,y0,z0])

    rot = np.zeros(data.size, dtype=data.dtype)
    n = data.size
    for name in data.dtype.names:
        rot[name][:n] = data[name]

    rot = add_col(rot, 'x0', np.array([x0]*data.size), dtype=float)
    rot = add_col(rot, 'y0', np.array([y0]*data.size), dtype=float)
    rot = add_col(rot, 'z0', np.array([z0]*data.size), dtype=float)

    for i in xrange(n):

        # We can skip the rotation for central objects
        if (data['most_massive'][i]==1):
            if verbose:
                print 'skipping object -- it is classified as a central galaxy'
            pos = positions.T[i]
            R = np.sqrt(sum(pos*pos))
            rot['rh'][i] = R
            continue

        if verbose:
            print i,

        # Cartesian galaxy position, relative to the halo centre
        pos = positions.T[i]

        # Generate a new pair of randomised position angles on a sphere of radius R
        phi_shifted = np.random.rand() * np.pi # range: [0, \pi]
        theta_shifted = (np.random.rand()-0.5) * 2 * np.pi # range: [-\pi, \pi]

        # Galaxy position in radial coordinates
        R = np.sqrt(sum(pos*pos))
        phi = np.arccos(pos[2]/R) * pos[0]/abs(pos[0])
        theta = np.arcsin(pos[1]/R/np.sin(phi))


        #import pdb ; pdb.set_trace()
        Rxyz = build_rotation_matrix()



       # g0 = np.array([ [np.dot(pos,rotated), -get_norm(np.cross(pos,rotated)), 0], [get_norm(np.cross(pos,rotated)), np.dot(pos,rotated),  0], [0,0,1]])
       # f0 = np.array([ pos, (rotated-np.dot(pos,rotated)*pos)/get_norm(rotated-np.dot(pos,rotated)*pos), np.cross(rotated,pos) ])
       # Rmat = f0*g0/f0

        # Work out the rotated galaxy position
        #xrot = R * np.sin(phi_shifted) * np.cos(theta_shifted)
        #yrot = R * np.sin(phi_shifted) * np.sin(theta_shifted)
        #zrot = R * np.cos(phi_shifted)

        #rotated = np.array([xrot, yrot, zrot])

        # Then do the 3D shape vector
        #Rxyz = np.array([[np.cos(dtheta)*np.cos(dpsi), -np.cos(dphi)*np.sin(dpsi) + np.sin(dphi)*np.sin(dtheta)*np.cos(dpsi), np.sin(dpsi)*np.sin(dphi) + np.sin(dtheta)*np.cos(dpsi)*np.cos(dphi) ],
        #    [np.sin(dpsi)*np.cos(dtheta), np.cos(dphi)*np.cos(dpsi) + np.sin(dphi)*np.sin(dtheta)*np.sin(dpsi), -np.sin(dphi)*np.cos(dpsi) + np.cos(dphi)*np.sin(dphi)*np.sin(dpsi)],
        #    [-np.sin(dtheta), np.sin(dphi)*np.cos(dtheta), np.cos(dtheta)*np.cos(dphi)]])

        rotated = np.dot(Rxyz,pos)
        xrot = rotated[0]
        yrot = rotated[1]
        zrot = rotated[2]

        if verbose:
         #   print 'ROTATING THROUGH 180 deg'
            print 'New position (x, y, z) : %3.3f, %3.3f %3.3f'%(xrot,yrot,zrot)

        a3d = np.array([data['a1'][i], data['a2'][i], data['a3'][i]])
        b3d = np.array([data['b1'][i], data['b2'][i], data['b3'][i]])
        c3d = np.array([data['c1'][i], data['c2'][i], data['c3'][i]])
        arot = np.dot(Rxyz,a3d)
        brot = np.dot(Rxyz,b3d)
        crot = np.dot(Rxyz,c3d)

        a3d_dm = np.array([data['a1_dm'][i], data['a2_dm'][i], data['a3_dm'][i]])
        b3d_dm = np.array([data['b1_dm'][i], data['b2_dm'][i], data['b3_dm'][i]])
        c3d_dm = np.array([data['c1_dm'][i], data['c2_dm'][i], data['c3_dm'][i]])
        arot_dm = np.dot(Rxyz,a3d_dm)
        brot_dm = np.dot(Rxyz,b3d_dm)
        crot_dm = np.dot(Rxyz,c3d_dm)

        #import pdb ; pdb.set_trace()

        #arot = #rotate_shape_vector(data[i], theta_shifted , phi_shifted, theta, phi, axis='a%d')
        #brot = #rotate_shape_vector(data[i], theta_shifted , phi_shifted, theta, phi, axis='b%d')
        #crot = #rotate_shape_vector(data[i], theta_shifted , phi_shifted, theta, phi, axis='c%d')

        #arot_dm = #rotate_shape_vector(data[i], theta_shifted , phi_shifted, theta, phi, axis='a%d_dm')
        #brot_dm = #rotate_shape_vector(data[i], theta_shifted , phi_shifted, theta, phi, axis='b%d_dm')
        #crot_dm = #rotate_shape_vector(data[i], theta_shifted , phi_shifted, theta, phi, axis='c%d_dm')
        
        rot['npart_baryon'][i] = data['npart_baryon'][i]
        rot['rh'][i] = R
        
        rot['x'][i] = copy.deepcopy(rotated[0])+x0
        rot['y'][i] = copy.deepcopy(rotated[1])+y0
        rot['z'][i] = copy.deepcopy(rotated[2])+z0

        # Need to shift the objects left outside the box by the rotation back
        # to the other side if the universe 
        if verbose: 
            print 'Rotation leaves %d objects outside the simulation box.'%(rot['x'][(rot['x']<0) | (rot['y']<0) | (rot['z']<0)].size)

        if (rot['x'][i]<0):
            rot['x'][i]+=100
        if (rot['y'][i]<0):
            #print 'Shifting in y %3.3f'%rot['y'][i]
            rot['y'][i]+=100
        if (rot['z'][i]<0):
            #print 'Shifting in z %3.3f'%rot['z'][i]
            rot['z'][i]+=100
        if (rot['x'][i]>100):
            rot['x'][i]-=100
        if (rot['y'][i]>100):
            rot['y'][i]-=100
        if (rot['z'][i]>100):
            rot['z'][i]-=100
       
        rot['c3'][i]=crot[2]
        rot['c2'][i]=crot[1]
        rot['c1'][i]=crot[0]
        rot['b3'][i]=brot[2]
        rot['b2'][i]=brot[1]
        rot['b1'][i]=brot[0]
        rot['a3'][i]=arot[2]
        rot['a2'][i]=arot[1]
        rot['a1'][i]=arot[0]

        rot['c3_dm'][i]=crot_dm[2]
        rot['c2_dm'][i]=crot_dm[1]
        rot['c1_dm'][i]=crot_dm[0]
        rot['b3_dm'][i]=brot_dm[2]
        rot['b2_dm'][i]=brot_dm[1]
        rot['b1_dm'][i]=brot_dm[0]
        rot['a3_dm'][i]=arot_dm[2]
        rot['a2_dm'][i]=arot_dm[1]
        rot['a1_dm'][i]=arot_dm[0]

        # Work out the projected 2D ellipticities
        a3d = np.array([[rot['a1'][i], rot['a2'][i], rot['a3'][i]]])
        b3d = np.array([[rot['b1'][i], rot['b2'][i], rot['b3'][i]]])
        c3d = np.array([[rot['c1'][i], rot['c2'][i], rot['c3'][i]]])
        q3d = np.array([np.sqrt(rot['lambda_a'][i]/rot['lambda_c'][i])])
        s3d = np.array([np.sqrt(rot['lambda_b'][i]/rot['lambda_c'][i])])
        e1,e2 = project_3d_shape(a3d, b3d, c3d, q3d, s3d)
        rot['e1'][i] = e1
        rot['e2'][i] = e2

        #import pdb ; pdb.set_trace()

        

        #if abs(np.sqrt(e1*e1+e2*e2))>0.75:
        #    import pdb ; pdb.set_trace()
        


        # Work out the projected 2D ellipticities
        a3d = np.array([[rot['a1_dm'][i], rot['a2_dm'][i], rot['a3_dm'][i]]])
        b3d = np.array([[rot['b1_dm'][i], rot['b2_dm'][i], rot['b3_dm'][i]]])
        c3d = np.array([[rot['c1_dm'][i], rot['c2_dm'][i], rot['c3_dm'][i]]])
        q3d = np.array([np.sqrt(rot['lambda_a_dm'][i]/rot['lambda_c_dm'][i])])
        s3d = np.array([np.sqrt(rot['lambda_b_dm'][i]/rot['lambda_c_dm'][i])])
        e1,e2 = project_3d_shape(a3d, b3d, c3d, q3d, s3d)
        rot['e1_dm'][i] = e1
        rot['e2_dm'][i] = e2

    return rot




def rotate_shape_vector(data, theta_shifted, phi_shifted, theta, phi, axis='a%d'):
    # Axis vector in Cartesian, then radial coordinates
    # assume it's normalised correctly in 3D
    avec = np.array([data[axis%1], data[axis%2], data[axis%3]])
    #phi_s = np.arccos(avec[2]) * avec[0]/abs(avec[0])
    theta_s = np.arctan2(avec[1],avec[0]) #np.arcsin(avec[1]/np.sin(phi_s))
    axy=np.sqrt(avec[0]*avec[0]+avec[1]*avec[1])
    phi_s = np.pi - np.arctan2(avec[2], axy)

    # Work out the correct new shape vector, preserving the relative orientation to the centre of mass
    dphi = phi_s - phi
    dtheta = theta_s - theta

    phi_s_shifted = phi_shifted + dphi
    theta_s_shifted = theta_shifted + dtheta

    axy0 = np.sin(phi_s_shifted)
    ax_shifted = axy0 * np.cos(theta_s_shifted)
    ay_shifted = axy0 * np.sin(theta_s_shifted)
    az_shifted = axy0 * np.tan(np.pi - phi_s_shifted)

    rot = np.array([ax_shifted, ay_shifted, az_shifted])

    return rot

def evaluate_ellipticity(cat):
    # Work out the ellipticities from the projected axis lengths for a catalogue of galaxies 
    phi = np.arctan2(cat['a2_2d'], cat['a1_2d'])
    q = np.sqrt(cat['lambda2_2d']/cat['lambda1_2d'])
    e = (q-1)/(q+1)
    e[np.invert(np.isfinite(e))] = 0.0

    return e, e*np.cos(2*phid), e*np.sin(2*phid)




def symmetrise_halo3(data, gids=None, verbose=True, g=None):
    if verbose:
        print 'Will apply one random rotation per subhalo'
    np.random.seed(9000)

    # Work out the centroid about which to rotate
    import weightedstats as ws
    N = data['npart_baryon']
    # Let's try something slightly different here
    # get the centroid position from the database rather than trying to recalculate it
    info = find_centre(g)

    #info = find_mean_position(g, data)
    x0, y0, z0 = info['x']/1000., info['y']/1000., info['z']/1000.

    cent = np.array([x0,y0,z0])
    positions = np.array([data['x']-x0, data['y']-y0, data['z']-z0])
    c = np.array([data['c3'], data['c2'], data['c1']])
    c_dm = np.array([data['c3_dm'], data['c2_dm'], data['c1_dm']])
    a = np.array([data['a3'], data['a2'], data['a1']])
    a_dm = np.array([data['a3_dm'], data['a2_dm'], data['a1_dm']])
    b = np.array([data['b3'], data['b2'], data['b1']])
    b_dm = np.array([data['b3_dm'], data['b2_dm'], data['b1_dm']])

    rot = np.zeros(data.size, dtype=data.dtype)
    n = data.size
    for name in data.dtype.names:
        rot[name][:n] = data[name]

    for i in xrange(n):

        # We can skip the rotation for central objects
        if (data['central'][i]==1):
            print 'skipping object -- it is classified as a central galaxy'
            continue

        if verbose:
            print i,
        # Define a random rotation about each axis
        alpha = 2 * np.pi * (np.random.rand() - 0.5) 
        beta = 2 * np.pi * (np.random.rand() - 0.5) 
        gamma = 2 * np.pi * (np.random.rand() - 0.5)

        if verbose:
            print alpha, beta, gamma

        pos = positions.T[i]


        # Do the symmetrisation as three independent rotations in the xy, xz and yz planes 
        rotated = copy.deepcopy(pos)
        Rxy = np.array([[np.cos(alpha), -np.sin(alpha)], [np.sin(alpha), np.cos(alpha)]])
        Rxz = np.array([[np.cos(beta), -np.sin(beta)], [np.sin(beta), np.cos(beta)]])
        Ryz = np.array([[np.cos(gamma), -np.sin(gamma)], [np.sin(gamma), np.cos(gamma)]])
        rotated[:2] = np.dot(Rxy,rotated[:2])
        rotated[[0,2]] = np.dot(Rxz,np.array([rotated[0],rotated[2]]))
        rotated[1:] = np.dot(Ryz,rotated[1:])


        # Now do the shapes in the same way
        arot = copy.deepcopy(a.T[i])
        arot[:2] = np.dot(Rxy,arot[:2])
        arot[[0,2]] = np.dot(Rxz,np.array([arot[0],arot[2]]))
        arot[1:] = np.dot(Ryz,arot[1:])

        brot = copy.deepcopy(b.T[i])
        brot[:2] = np.dot(Rxy,brot[:2])
        brot[[0,2]] = np.dot(Rxz,np.array([brot[0],brot[2]]))
        brot[1:] = np.dot(Ryz,brot[1:])

        crot = copy.deepcopy(c.T[i])
        crot[:2] = np.dot(Rxy,crot[:2])
        crot[[0,2]] = np.dot(Rxz,np.array([crot[0],crot[2]]))
        crot[1:] = np.dot(Ryz,crot[1:])

        arot_dm = copy.deepcopy(a_dm.T[i])
        arot_dm[:2] = np.dot(Rxy,arot_dm[:2])
        arot_dm[[0,2]] = np.dot(Rxz,np.array([arot_dm[0],arot_dm[2]]))
        arot_dm[1:] = np.dot(Ryz,arot_dm[1:])

        brot_dm = copy.deepcopy(b_dm.T[i])
        brot_dm[:2] = np.dot(Rxy,brot_dm[:2])
        brot_dm[[0,2]] = np.dot(Rxz,np.array([brot_dm[0],brot_dm[2]]))
        brot_dm[1:] = np.dot(Ryz,brot_dm[1:])

        crot_dm = copy.deepcopy(c_dm.T[i])
        crot_dm[:2] = np.dot(Rxy,crot_dm[:2])
        crot_dm[[0,2]] = np.dot(Rxz,np.array([crot_dm[0],crot_dm[2]]))
        crot_dm[1:] = np.dot(Ryz,crot_dm[1:])

        rot['npart_baryon'][i] = data['npart_baryon'][i]

        
        rot['x'][i] = copy.deepcopy(rotated[0])+x0
        rot['y'][i] = copy.deepcopy(rotated[1])+y0
        rot['z'][i] = copy.deepcopy(rotated[2])+z0

        #import pdb ; pdb.set_trace()

        rot['c3'][i]=crot[0]
        rot['c2'][i]=crot[1]
        rot['c1'][i]=crot[2]
        rot['b3'][i]=brot[0]
        rot['b2'][i]=brot[1]
        rot['b1'][i]=brot[2]
        rot['a3'][i]=arot[0]
        rot['a2'][i]=arot[1]
        rot['a1'][i]=arot[2]

        rot['c3_dm'][i]=crot_dm[0]
        rot['c2_dm'][i]=crot_dm[1]
        rot['c1_dm'][i]=crot_dm[2]
        rot['b3_dm'][i]=brot_dm[0]
        rot['b2_dm'][i]=brot_dm[1]
        rot['b1_dm'][i]=brot_dm[2]
        rot['a3_dm'][i]=arot_dm[0]
        rot['a2_dm'][i]=arot_dm[1]
        rot['a1_dm'][i]=arot_dm[2]

        # Also rotate the ellipticities
        phim = np.arctan2(rot['a2'][i], rot['a1'][i])
        em = np.sqrt(rot['e1'][i]*rot['e1'][i]+rot['e2'][i]*rot['e2'][i])
        rot['e1'][i] = em * np.cos(2*phim)
        rot['e2'][i] = em * np.sin(2*phim)

        phid = np.arctan2(rot['a2_dm'][i], rot['a1_dm'][i])
        ed = np.sqrt(rot['e1_dm'][i]*rot['e1_dm'][i]+rot['e2_dm'][i]*rot['e2_dm'][i])
        rot['e1_dm'][i] = ed * np.cos(2*phid)
        rot['e2_dm'][i] = ed * np.sin(2*phid)


    return rot




def symmetrise_catalogue2(data=None,mask=None,filename='/home/ssamurof/massive_black_ii/subhalo_cat-nthreshold5.fits', savedir=None, replace=False, verbose=True, nmin=2, rank=0, size=1):
    root_folder='/physics/yfeng1/mb2'
    snapshot='085'

    if data is None:
        data = fi.FITS(filename)['baryons'].read()
        dm = fi.FITS(filename)['dm'].read()

        #mask = (data['npart']>0) & (dm['npart']>1000) & np.isfinite(data['x']) & np.isfinite(data['y'])  & np.isfinite(data['z'])

    if savedir is None:
        savedir = filename.replace('/subhalo_cat','/cats/symmetrised/subhalo_cat').replace('.fits', '-symmetrised-%d-masswtdmedian-stellarmasscut.fits'%rank)
    print savedir
    if savedir==filename:
        print "WARNING: target and input file paths are the same"
        import pdb ; pdb.set_trace()

    snap = SnapDir(snapshot, root_folder)
    h = snap.readsubhalo()

    group_ids = np.array(h['groupid'])[mask]
    groups = np.unique(group_ids)

    dt = np.dtype([('subhalo_id', int),('group_id', int), ('x', '>f8'), ('y', '>f8'), ('z', '>f8'), ('npart', '>f8'), ('lambda1', '>f8'), ('lambda2', '>f8'), ('lambda3', '>f8'), ('a1', '>f8'), ('a2', '>f8'), ('a3', '>f8'), ('b1', '>f8'), ('b2', '>f8'), ('b3', '>f8'), ('c1', '>f8'), ('c2', '>f8'), ('c3', '>f8')])

    # Copy over the old data
    outdat = np.zeros(len(data), dtype=dt)
    for name in dt.names:
        if (name=='subhalo_id') or (name=='group_id'):
            continue
        else:
            outdat[name] = data[name]

    ident = np.arange(0, len(data), 1)
    i0=0
    for i, g in enumerate(groups):

        if i%size!=rank:
            continue
        
        select = (group_ids==g)

        ngal = len(data[select])
        outdat['group_id'][select] = group_ids[select]
        if ngal<nmin:
            print 'Skipping halo %d as it contains <%d galaxies'%(g,nmin)
            outdat['subhalo_id'][select] = ident[select]
            continue

        #import pdb ; pdb.set_trace()

        symmetrised_halo = symmetrise_halo2(data[select], gids=group_ids[select], verbose=False, g=g)

        if verbose:
            print g, ngal

        for name in data.dtype.names:
            outdat[name][select] = symmetrised_halo[name]

        outdat['subhalo_id'][select] = ident[select]

        i0+=ngal

    outfits = fi.FITS(savedir, 'rw')
    outfits.write(outdat)
    outfits[-1].write_key('EXTNAME','baryons')

    print 'Done'
    return 0




def symmetrise_catalogue(data, mask=None, nrot=1000, savedir='/home/ssamurof/massive_black_ii/cats/symmetrised/', replace=False):
    root_folder='/physics/yfeng1/mb2'
    snapshot='085'

    savedir+='nrot%d'%nrot

    files = glob.glob('%s/*fits'%savedir)
    groups_done = []
    if not replace and (len(files)>0):
        for f in files:
            fits = fi.FITS(f) 
            print f, len(fits)
            for hdu in fits[1:]:
                groups_done.append(hdu.read_header()['EXTNAME'])
        groups_done = [g.replace(' ', '') for g in groups_done] 

    if not os.path.exists(savedir):
        os.system('mkdir -p %s'%savedir)

    snap = SnapDir(snapshot, root_folder)
    h = snap.readsubhalo()

    if mask is None:
        mask = np.one(data.size).astype(bool)

    group_ids = np.array(h['groupid'])[mask]
    groups = np.unique(group_ids)

    outfits = fi.FITS('%s/symmetrised_subhalo_cat-halo%d.fits'%(savedir,groups[0]),'rw')

    for i, g in enumerate(groups):
        
        # Check against the list of halo ids that have already been processed
        if (not replace) and ('halo%d'%g in groups_done):
            continue
        else:
            print i, g

        select = (group_ids==g)

        symmetrised_halo = symmetrise_halo(data[mask][select], gids=group_ids[select], nrot=nrot, verbose=False)

        if (i%5000==0) and (i!=0):
            outfits.close()
            outfits = fi.FITS('%s/symmetrised_subhalo_cat-halo%d.fits'%(savedir,g),'rw')

        outfits.write(symmetrised_halo)
        outfits[-1].write_key('EXTNAME','halo%d'%g)

    print 'Done'
    return 0



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

def find_mean_position(g, data):
    info = np.zeros(1, dtype=[('x',float), ('y',float), ('z',float)])

    for name in ['x','y','z']:
        info[name] = np.mean(data[name][data['halo_id']==g])

    return info
    

def symmetrise_halo2(data, gids=None, verbose=True, g=None):
    if verbose:
        print 'Will apply one random rotation per subhalo'
    np.random.seed(9000)

    # Work out the centroid about which to rotate
    import weightedstats as ws
    N = data['npart']
    # Let's try something slightly different here
    # get the centroid position from the database rather than trying to recalculate it
    info = find_centre(g) 
    x0, y0, z0 = info['x'], info['y'], info['z']

#    xrand = np.random.choice(data['x'])
#    yrand = np.random.choice(data['y'])
#    zrand = np.random.choice(data['z'])
#    sane = N>0 # stellar mass cut 
#    #(abs(data['x']-xrand)<0.1e5) & (abs(data['y']-yrand)<0.1e5) & (abs(data['z']-zrand)<0.1e5) 
#    x0 = ws.numpy_weighted_median(data['x'][np.isfinite(data['x']) & sane], weights=N[sane & np.isfinite(data['x'])])
#    y0 = ws.numpy_weighted_median(data['y'][np.isfinite(data['y']) & sane], weights=N[sane & np.isfinite(data['y'])])
#    z0 =  ws.numpy_weighted_median(data['z'][np.isfinite(data['z']) & sane], weights=N[sane & np.isfinite(data['z'])])
    cent = np.array([x0,y0,z0])
    positions = np.array([data['x']-x0, data['y']-y0, data['z']-z0])
    c = np.array([data['c3'], data['c2'], data['c1']])
    a = np.array([data['a3'], data['a2'], data['a1']])
    b = np.array([data['b3'], data['b2'], data['b1']])

    rot = np.zeros(data.size, dtype=data.dtype)
    n = data.size
    for name in data.dtype.names:
        rot[name][:n] = data[name]

    for i in xrange(n):
        if verbose:
            print i,
        # Define a random rotation about each axis
        alpha = 2 * np.pi * (np.random.rand() - 0.5) 
        beta = 2 * np.pi * (np.random.rand() - 0.5) 
        gamma = 2 * np.pi * (np.random.rand() - 0.5)

        if verbose:
            print alpha, beta, gamma
        #alpha,beta,gamma=np.pi,0,0

        # Rotate the position vector
        #Rx = np.array([[1,0,0], [0, np.cos(alpha), -np.sin(alpha)], [0, np.sin(alpha), np.cos(alpha)]])
        #rotated = np.dot(Rx,positions.T[i])
        #Ry = np.array([[np.cos(beta), 0, np.sin(beta)],[0, 1, 0],[-np.sin(beta), 0, np.cos(beta)]])
        #rotated = np.dot(Ry,rotated)
        #Rz = np.array([[np.cos(gamma), -np.sin(gamma), 0], [np.sin(gamma), np.cos(gamma), 0], [0,0,1]])
        #rotated = np.dot(Rz,rotated)

        pos = positions.T[i]


        # Do the symmetrisation as three independent rotations in the xy, xz and yz planes 
        rotated = copy.deepcopy(pos)
        Rxy = np.array([[np.cos(alpha), -np.sin(alpha)], [np.sin(alpha), np.cos(alpha)]])
        Rxz = np.array([[np.cos(beta), -np.sin(beta)], [np.sin(beta), np.cos(beta)]])
        Ryz = np.array([[np.cos(gamma), -np.sin(gamma)], [np.sin(gamma), np.cos(gamma)]])
        rotated[:2] = np.dot(Rxy,rotated[:2])
        rotated[[0,2]] = np.dot(Rxz,np.array([rotated[0],rotated[2]]))
        rotated[1:] = np.dot(Ryz,rotated[1:])


        # Now do the shapes in the same way
        arot = copy.deepcopy(a.T[i])
        arot[:2] = np.dot(Rxy,arot[:2])
        arot[[0,2]] = np.dot(Rxz,np.array([arot[0],arot[2]]))
        arot[1:] = np.dot(Ryz,arot[1:])

        brot = copy.deepcopy(b.T[i])
        brot[:2] = np.dot(Rxy,brot[:2])
        brot[[0,2]] = np.dot(Rxz,np.array([brot[0],brot[2]]))
        brot[1:] = np.dot(Ryz,brot[1:])

        crot = copy.deepcopy(c.T[i])
        crot[:2] = np.dot(Rxy,crot[:2])
        crot[[0,2]] = np.dot(Rxz,np.array([crot[0],crot[2]]))
        crot[1:] = np.dot(Ryz,crot[1:])

        rot['npart'][i] = data['npart'][i]

        
        rot['x'][i] = copy.deepcopy(rotated[0])+x0
        rot['y'][i] = copy.deepcopy(rotated[1])+y0
        rot['z'][i] = copy.deepcopy(rotated[2])+z0

        rot['c3'][i]=crot[0]
        rot['c2'][i]=crot[1]
        rot['c1'][i]=crot[2]
        rot['b3'][i]=brot[0]
        rot['b2'][i]=brot[1]
        rot['b1'][i]=brot[2]
        rot['a3'][i]=arot[0]
        rot['a2'][i]=arot[1]
        rot['a1'][i]=arot[2]

    #import pdb ; pdb.set_trace()

    return rot


def symmetrise_halo(data, nrot=100, gids=None, verbose=True):
    if verbose:
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
        if verbose:
            print i,
        # Define a random rotation about each axis
        alpha = 2 * np.pi * (np.random.rand() - 0.5) 
        beta = 2 * np.pi * (np.random.rand() - 0.5) 
        gamma = 2 * np.pi * (np.random.rand() - 0.5)

        if verbose:
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


def construct_random_cat(data, mask=None, format='halotools', f=1, fixed_bounds=True):
    """Construct a catalogue of random points drawn from the same volume as the data provided."""

    if f!=1:
        print 'Random point multiplier: %3.3f'%f
    if mask is None:
        mask = np.ones(data.size).astype(bool)
    if not fixed_bounds:
        rx = (np.random.random(size=int(data['x'][mask].size*f)) - 0.5) * (data['x'][mask].max()-data['x'][mask].min()) + data['x'][mask].mean()
        ry = (np.random.random(size=int(data['x'][mask].size*f)) - 0.5) * (data['y'][mask].max()-data['y'][mask].min()) + data['y'][mask].mean()
        rz = (np.random.random(size=int(data['x'][mask].size*f)) - 0.5) * (data['z'][mask].max()-data['z'][mask].min()) + data['z'][mask].mean()
    else:
        print 'Enforcing box edges'
        rx = np.random.random(size=int(data['x'][mask].size*f)) * 100000.
        ry = np.random.random(size=int(data['x'][mask].size*f)) * 100000.
        rz = np.random.random(size=int(data['x'][mask].size*f)) * 100000.

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

def gather_mpi_output(filestring, hdu=-1, name='subhalo_id', save=False):
    files = glob.glob(filestring)
    print 'Found %d files'%len(files)
    dat = fi.FITS(files[0])[hdu].read()
    for f in files:
        print f,
        info = fi.FITS(f)[hdu].read()
        mask = (info[name]!=0)
        print '%d objects'%len(info[name][mask])
        dat[mask] = info[mask]

    if save:
        print 'Saving results as single FITS file.'
        out = fi.FITS('bundled_fits.fits', 'rw')
        out.write(dat)
        if type(hdu) is str:
            out[-1].write_key('EXTNAME', hdu)
        out.close()
    return dat

def export_treecorr_output(filename,corr,errors):
    R = np.exp(corr.logr)

    xi = corr.xi

    out = np.vstack((R,xi,corr.weight, errors))
    print 'Saving %s'%filename
    np.savetxt(filename, out.T)

# This script module contains a number of utility functions
import numpy as np
import treecorr as tc

def build_rotation_matrix():
    """Generates a random rotation matrix, 
    which transforms a given 3D position vector to a random position on the sphere.
    Usage : rotated = R.unrotated"""
    u1 = np.random.rand()
    u2 = np.random.rand()

    # Use two random values from a uniform distribution to generate an axis about which to rotate
    theta = np.arccos(1*u1-1) 
    phi = 2 * np.pi * u2

    # Work out a Cartesian unit vector defined by the rotation axis 
    x = np.sin(theta)*np.cos(phi)
    y = np.sin(theta)*np.sin(phi)
    z = np.cos(theta)
    vec= np.array([x,y,z])
    vec/=get_norm(vec)

    # Now generate a random rotation angle about that axis
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




def get_norm(vec):
    return np.sqrt(sum([v*v for v in vec])) 

def project_3d_shape(a3d, b3d, c3d, q3d, s3d):
    """
    This function projects 3D ellipsoidal shapes onto a 2D ellipse and returns
    the 2 cartesian components of the ellipticity.

    See Piras2017:1707.06559 section 3.1
    and Joachimi2013:1203.6833
    """

    s = np.stack([a3d, b3d, c3d])
    w = np.stack([np.ones_like(q3d), q3d, s3d])

    k = np.sum(s[:,:,0:2]*np.expand_dims(s[:,:,2], axis=-1) / np.expand_dims(w[:,:]**2, axis=-1), axis=0)
    a2 =np.sum(s[:,:,2]**2/w[:,:]**2, axis=0)
    Winv = np.sum(np.einsum('ijk,ijl->ijkl', s[:,:,0:2], s[:,:,0:2]) /
               np.expand_dims(np.expand_dims(w[:,:]**2,-1),-1), axis=0) - np.einsum('ij,ik->ijk', k,k)/np.expand_dims(np.expand_dims(a2,-1),-1)
    W = np.linalg.inv(Winv)
    d = np.sqrt(np.linalg.det(W))
    e1 = (W[:,0,0] - W[:,1,1])/( W[:,0,0] + W[:,1,1] + 2*d)
    e2 = 2 * W[:,0,1]/( W[:,0,0] + W[:,1,1] + 2*d)
    return e1, e2

def project_to_sky(data):
    r = np.sqrt(data['x']*data['x'] + data['y']*data['y'] + data['z']*data['z'])
    ra = np.arccos(data['z']/r)
    dec = np.arccos(data['x']/r/np.sin(ra))

    ra = ra * 180./np.pi
    dec = dec * 180./np.pi

    fig = plt.figure(0)
    footprint_sub(ra, dec, 10, 10, 2048, fig, cmap='Purples')

    return fig

