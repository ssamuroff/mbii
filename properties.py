import numpy as np
from readsubhalo import *
import fitsio as fi

def compute_inertia_tensors(snap, reduced=False, inclusion_threshold=5):
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
    
    eigvectors = np.zeros((2, len(h), 3, 3))
    eigvalues  = np.zeros((2, len(h), 3))
    centroids  = np.zeros((2, len(h), 3))
    length  = np.zeros((2, len(h)))
 
    # Will compute for each halo the inertia tensor
    for i in range(len(h)):
        if i%100 ==0:
            print "Done %d samples"%i
        # Reject subhalos with less than some threshold occupation number
        if len(x[i]) < inclusion_threshold:
            pass  
        else: 
            weights = np.ones(len(x[i]))
            normFactor = np.double(len(x[i].T[0]))

            x0 = x[i].T[0] - np.dot(x[i].T[0],weights)/normFactor
            x1 = x[i].T[1] - np.dot(x[i].T[1],weights)/normFactor
            x2 = x[i].T[2] - np.dot(x[i].T[2],weights)/normFactor

            # Beware the ambiguous notation here - 
            # These weights define the difference between reduced and standard inertia tensors.
            # i.e. downweighting particles at the fringes of the subhalo
            # They have nothing to do with whether the mass weighting is applied.
            if reduced:
                wt = (x0*x0 + x1*x1 + x2*x2)
                select = (wt!=0)
            else:
                wt = np.ones(x0.size)
                select = wt.astype(bool)

            normFactor = np.double(len(x[i].T[0][select]))

            tens = np.zeros((3,3))
            tens[0, 0] = np.dot(weights[select] * x0[select], x0[select]/wt[select])/normFactor
            tens[1, 1] = np.dot(weights[select] * x1[select], x1[select]/wt[select])/normFactor
            tens[2, 2] = np.dot(weights[select] * x2[select], x2[select]/wt[select])/normFactor
            tens[1, 0] = np.dot(weights[select] * x1[select], x0[select]/wt[select])/normFactor
            tens[2, 0] = np.dot(weights[select] * x2[select], x0[select]/wt[select])/normFactor
            tens[2, 1] = np.dot(weights[select] * x2[select], x1[select]/wt[select])/normFactor
            tens[0, 2] = tens[2, 0]
            tens[0, 1] = tens[1, 0]
            tens[1, 2] = tens[2, 1]

            # Evaluate the mass-weighted centroid along each axis
            X = np.trapz(m[i]*x[i].T[0], x[i].T[0])/np.trapz(m[i], x[i].T[0])
            Y = np.trapz(m[i]*x[i].T[1], x[i].T[1])/np.trapz(m[i], x[i].T[1])
            Z = np.trapz(m[i]*x[i].T[2], x[i].T[2])/np.trapz(m[i], x[i].T[2])

            # Compute the eigenvalues of the halos and store the outputs
            w, v = np.linalg.eigh(tens)
            eigvalues[0,i] = w
            eigvectors[0,i] = v
            length[0,i] = normFactor
            centroids[0,i] = np.array([X,Y,Z])

        if (len(xb[i]) < inclusion_threshold):
            pass
        else:

            weights = np.ones(len(xb[i]))
            normFactor = np.double(len(xb[i].T[0]))

            x0 = xb[i].T[0] - np.dot(xb[i].T[0],weights)/normFactor
            x1 = xb[i].T[1] - np.dot(xb[i].T[1],weights)/normFactor
            x2 = xb[i].T[2] - np.dot(xb[i].T[2],weights)/normFactor

            if reduced:
                wt = (x0*x0 + x1*x1 + x2*x2)
                select = (wt!=0)
            else:
                wt = np.ones(x0.size)
                select = wt.astype(bool)

            normFactor = np.double(len(xb[i].T[0][select]))
            
            tens = np.zeros((3,3))
            tens[0, 0] = np.dot(weights[select] * x0[select], x0[select]/wt[select])/normFactor
            tens[1, 1] = np.dot(weights[select] * x1[select], x1[select]/wt[select])/normFactor
            tens[2, 2] = np.dot(weights[select] * x2[select], x2[select]/wt[select])/normFactor
            tens[1, 0] = np.dot(weights[select] * x1[select], x0[select]/wt[select])/normFactor
            tens[2, 0] = np.dot(weights[select] * x2[select], x0[select]/wt[select])/normFactor
            tens[2, 1] = np.dot(weights[select] * x2[select], x1[select]/wt[select])/normFactor
            tens[0, 2] = tens[2, 0]
            tens[0, 1] = tens[1, 0]
            tens[1, 2] = tens[2, 1]

            # Evaluate the mass-weighted centroid along each axis
            X = np.trapz(mb[i]*xb[i].T[0], xb[i].T[0])/np.trapz(mb[i], xb[i].T[0])
            Y = np.trapz(mb[i]*xb[i].T[1], xb[i].T[1])/np.trapz(mb[i], xb[i].T[1])
            Z = np.trapz(mb[i]*xb[i].T[2], xb[i].T[2])/np.trapz(mb[i], xb[i].T[2])

            # Compute the eigenvalues of the halos
            try:
                w, v = np.linalg.eigh(tens)
            except:
                import pdb ; pdb.set_trace()
            eigvalues[1,i] = w
            eigvectors[1,i] = v
            length[1,i] = normFactor
            centroids[1,i] = np.array([X,Y,Z])



    print "Saving output"
    if reduced:
        out = fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat_reduced-nthreshold%d.fits'%(inclusion_threshold),'rw')
    else:
        out = fi.FITS('/home/ssamurof/massive_black_ii/subhalo_cat-nthreshold%d.fits'%(inclusion_threshold),'rw')

    dat=np.zeros(len(h), dtype=[('x', float), ('y', float), ('z', float), ('npart',float), ('lambda1', float), ('lambda2', float), ('lambda3', float), ('a1', float), ('a2', float), ('a3', float), ('b1', float), ('b2', float), ('b3', float), ('c1', float), ('c2', float), ('c3', float)])
    dat['lambda1'] = eigvalues[0].T[0]
    dat['lambda2'] = eigvalues[0].T[1]
    dat['lambda3'] = eigvalues[0].T[2]
    dat['a1'] = eigvectors[0].T[0,0]
    dat['a2'] = eigvectors[0].T[0,1]
    dat['a3'] = eigvectors[0].T[0,2]
    dat['b1'] = eigvectors[0].T[1,0]
    dat['b2'] = eigvectors[0].T[1,1]
    dat['b3'] = eigvectors[0].T[1,2]
    dat['c1'] = eigvectors[0].T[2,0]
    dat['c2'] = eigvectors[0].T[2,1]
    dat['c3'] = eigvectors[0].T[2,2]
    dat['x'] = centroids[0].T[0]
    dat['y'] = centroids[0].T[1]
    dat['z'] = centroids[0].T[2]
    dat['npart'] = length[0]

    out.write(dat)
    out[-1].write_key('EXTNAME', 'dm')

    dat2=np.zeros(len(h), dtype=[('x', float), ('y', float), ('z', float), ('npart',float), ('lambda1', float), ('lambda2', float), ('lambda3', float), ('a1', float), ('a2', float), ('a3', float), ('b1', float), ('b2', float), ('b3', float), ('c1', float), ('c2', float), ('c3', float)])
    dat2['lambda1'] = eigvalues[1].T[0]      
    dat2['lambda2'] = eigvalues[1].T[1]
    dat2['lambda3'] = eigvalues[1].T[2]
    dat2['a1'] = eigvectors[1].T[0,0]
    dat2['a2'] = eigvectors[1].T[0,1]
    dat2['a3'] = eigvectors[1].T[0,2]
    dat2['b1'] = eigvectors[1].T[1,0]
    dat2['b2'] = eigvectors[1].T[1,1]
    dat2['b3'] = eigvectors[1].T[1,2]
    dat2['c1'] = eigvectors[1].T[2,0]
    dat2['c2'] = eigvectors[1].T[2,1]
    dat2['c3'] = eigvectors[1].T[2,2]
    dat2['x'] = centroids[1].T[0]
    dat2['y'] = centroids[1].T[1]
    dat2['z'] = centroids[1].T[2]
    dat2['npart'] = length[1]

    out.write(dat2)
    out[-1].write_key('EXTNAME', 'baryons')
    out.close()
import numpy as np
from readsubhalo import *
import fitsio as fi

def compute_spin(snap, inclusion_threshold=1, component='baryons', nsubhalo=-1):
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


