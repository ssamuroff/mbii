import numpy as np
from readsubhalo import *
import fitsio as fi

def compute_inertia_tensors_projected(snap, reduced=False, inclusion_threshold=1):
    """ Compute the intertia tensors for all subhalos for
    the dark matter and stellar components and saves the output for
    later integration within the database
    """

    print 'Using reduced (distance weighted) tensors', 'yes'*int(reduced), 'no'*int(np.invert(reduced) )

    # Read the subhalo information
    h = snap.readsubhalo()
    # Load the positions and masses of the constituent particles
    print 'Loading dark matter particles'
    x = snap.load(1, 'pos', h)
    m = snap.load(1, 'mass', h)
    print 'Loading baryon particles'
    xb = snap.load(4, 'pos', h)
    mb = snap.load(4, 'mass', h)
    
    eigvectors = np.zeros((2, len(h), 2, 2))
    eigvalues  = np.zeros((2, len(h), 2))
    centroids  = np.zeros((2, len(h), 3))
    spher_pos  = np.zeros((2, len(h), 3))
    length  = np.zeros((2, len(h)))
 
    # Will compute for each halo the inertia tensor
    for i in range(len(h)):
        if i%100 ==0:
            print "Done %d samples"%i
        # Reject subhalos with less than some threshold occupation number
        if len(x[i]) < inclusion_threshold:
            pass  
        else: 
            # Decide whether to weight by the particle mass
            weights = np.ones(len(x[i]))
            normFactor = np.double(len(x[i].T[0]))

            x0 = x[i].T[0] - np.dot(x[i].T[0],weights)/normFactor
            x1 = x[i].T[1] - np.dot(x[i].T[1],weights)/normFactor

            # Beware the ambiguous notation here - 
            # These weights define the difference between reduced and standard inertia tensors.
            # i.e. downweighting particles at the fringes of the subhalo
            # They have nothing to do with whether the mass weighting is applied.
            if reduced:
                wt = (x0*x0 + x1*x1)
                select = (wt!=0)
            else:
                wt = np.ones(x0.size)
                select = wt.astype(bool)

            normFactor = np.double(len(x[i].T[0][select]))

            tens = np.zeros((2,2))
            tens[0, 0] = np.dot(weights[select] * x0[select], x0[select]/wt[select])/normFactor
            tens[1, 1] = np.dot(weights[select] * x1[select], x1[select]/wt[select])/normFactor        
            tens[1, 0] = np.dot(weights[select] * x1[select], x0[select]/wt[select])/normFactor
            tens[0, 1] = tens[1, 0]

            # Evaluate the mass-weighted centroid along each axis
            X = np.trapz(m[i]*x[i].T[0], x[i].T[0])/np.trapz(m[i], x[i].T[0])
            Y = np.trapz(m[i]*x[i].T[1], x[i].T[1])/np.trapz(m[i], x[i].T[1])
            Z = np.trapz(m[i]*x[i].T[2], x[i].T[2])/np.trapz(m[i], x[i].T[2])

            phi = np.arctan2(Y,X)
            r = np.sqrt(X*X + Y*Y + Z*Z)
            theta = np.arccos(Z/r)

            # Compute the eigenvalues of the halos and store the outputs
            w, v = np.linalg.eigh(tens)
            eigvalues[0,i] = w
            eigvectors[0,i] = v
            spher_pos[0,i] = np.array([r,theta,phi])
            centroids[0,i] = np.array([X,Y,Z])
            length[0,i] = len(x[i])

        if (len(xb[i]) < inclusion_threshold):
            pass
        else:
            weights = np.ones(len(xb[i]))
            normFactor = np.double(len(xb[i].T[0]))

            x0 = xb[i].T[0] - np.dot(xb[i].T[0],weights)/normFactor
            x1 = xb[i].T[1] - np.dot(xb[i].T[1],weights)/normFactor
            if reduced:
                wt = (x0*x0 + x1*x1)
                select = (wt!=0)
            else:
                wt = np.ones(x0.size)
                select = wt.astype(bool)

            normFactor = np.double(len(xb[i].T[0][select]))

            tens = np.zeros((2,2))
            tens[0, 0] = np.dot(weights[select] * x0[select], x0[select]/wt[select])/normFactor
            tens[1, 1] = np.dot(weights[select] * x1[select], x1[select]/wt[select])/normFactor
            tens[1, 0] = np.dot(weights[select] * x1[select], x0[select]/wt[select])/normFactor
            tens[0, 1] = tens[1, 0]

            # Evaluate the mass-weighted centroid along each axis
            X = np.trapz(mb[i]*xb[i].T[0], xb[i].T[0])/np.trapz(mb[i], xb[i].T[0])
            Y = np.trapz(mb[i]*xb[i].T[1], xb[i].T[1])/np.trapz(mb[i], xb[i].T[1])
            Z = np.trapz(mb[i]*xb[i].T[2], xb[i].T[2])/np.trapz(mb[i], xb[i].T[2])

            phi = np.arctan2(Y,X)
            r = np.sqrt(X*X + Y*Y + Z*Z)
            theta = np.arccos(Z/r)

            # Compute the eigenvalues of the halos
            w, v = np.linalg.eigh(tens)
            eigvalues[1,i] = w
            eigvectors[1,i] = v
            spher_pos[1,i] = np.array([r,theta,phi])
            centroids[1,i] = np.array([X,Y,Z])
            length[1,i] = len(xb[i])


    print "Saving output"
    if reduced:
        out = fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat_reduced_projected-nthreshold%d.fits'%inclusion_threshold,'rw')
    else:
        out = fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat_projected-nthreshold%d.fits'%inclusion_threshold,'rw')

    dat=np.zeros(len(h), dtype=[('x', float), ('y', float), ('z', float), ('npart',float), ('lambda1', float), ('lambda2', float), ('a1', float), ('a2', float), ('b1', float), ('b2', float)])
    dat['lambda1'] = eigvalues[0].T[0]
    dat['lambda2'] = eigvalues[0].T[1]
    dat['a1'] = eigvectors[0].T[0,0]
    dat['a2'] = eigvectors[0].T[0,1]
    dat['b1'] = eigvectors[0].T[1,0]
    dat['b2'] = eigvectors[0].T[1,1]
    dat['x'] = centroids[0].T[0]
    dat['y'] = centroids[0].T[1]
    dat['z'] = centroids[0].T[2]
    dat['npart'] = length[0]

    out.write(dat)
    out[-1].write_key('EXTNAME', 'dm')

    dat2=np.zeros(len(h), dtype=[('x', float), ('y', float), ('z', float), ('npart',float), ('lambda1', float), ('lambda2', float), ('a1', float), ('a2', float), ('b1', float), ('b2', float)])
    dat2['lambda1'] = eigvalues[1].T[0]      
    dat2['lambda2'] = eigvalues[1].T[1]
    dat2['a1'] = eigvectors[1].T[0,0]
    dat2['a2'] = eigvectors[1].T[0,1]
    dat2['b1'] = eigvectors[1].T[1,0]
    dat2['b2'] = eigvectors[1].T[1,1]
    dat2['x'] = centroids[1].T[0]
    dat2['y'] = centroids[1].T[1]
    dat2['z'] = centroids[1].T[2]
    dat2['npart'] = length[1]

    out.write(dat2)
    out[-1].write_key('EXTNAME', 'baryons')
    out.close()