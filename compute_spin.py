import numpy as np
from readsubhalo import *
import fitsio as fi

def compute_inertia_tensors(snap):
    """ Compute the intertia tensors for all subhalos for
    the dark matter and stellar components and saves the output for
    later integration within the database
    """

    import pdb ; pdb.set_trace()

    # Read the subhalo information
    h = snap.readsubhalo()
    # Load the positions and masses of the constituent particles
    x = snap.load(1, 'pos', h)
    m = snap.load(1, 'mass', h)
    v = snap.load(1, 'vel', h)
    xb = snap.load(4, 'pos', h)
    mb = snap.load(4, 'mass', h)
    vb = snap.load(4, 'vel', h)
    
    eigvectors = np.zeros((2, len(h), 3, 3))
    eigvalues  = np.zeros((2, len(h), 3))
    centroids  = np.zeros((2, len(h), 3))
 
    # Will compute for each halo the inertia tensor
    for i in range(len(h)):
	if i%100 ==0:
	    print "Done %d samples"%i

	if len(x[i]) == 0:
	    pass # print "Halo %d is empty of dark matter"%i
	else:
	    weights = np.ones(len(x[i]))
	    normFactor = np.double(len(x[i].T[0]))
	    x0 = x[i].T[0] - np.dot(x[i].T[0],weights)/normFactor
	    x1 = x[i].T[1] - np.dot(x[i].T[1],weights)/normFactor
	    x2 = x[i].T[2] - np.dot(x[i].T[2],weights)/normFactor
	    tens = np.zeros((3,3))
	    tens[0, 0] = np.dot(x0,x0)/normFactor
	    tens[1, 1] = np.dot(x1,x1)/normFactor        
	    tens[2, 2] = np.dot(x2,x2)/normFactor
	    tens[1, 0] = np.dot(x1,x0)/normFactor
	    tens[2, 0] = np.dot(x2,x0)/normFactor
	    tens[2, 1] = np.dot(x2,x1)/normFactor
	    tens[0, 2] = tens[2, 0]
	    tens[0, 1] = tens[1, 0]
	    tens[1, 2] = tens[2, 1]
	    # Evaluate the mass-weighted centroid along each axis
	    X = np.trapz(m[i]*x[i].T[0], x[i].T[0])/np.trapz(m[i], x[i].T[0])
	    Y = np.trapz(m[i]*x[i].T[1], x[i].T[1])/np.trapz(m[i], x[i].T[1])
	    Z = np.trapz(m[i]*x[i].T[2], x[i].T[2])/np.trapz(m[i], x[i].T[2])
	    import pdb ; pdb.set_trace()
	    # Compute the eigenvalues of the halos
	    w, v = np.linalg.eigh(tens)
	    eigvalues[0,i] = w
	    eigvectors[0,i] = v
	    centroids[0,i] = np.array([X,Y,Z])
	if len(xb[i]) == 0:
	    pass #print "Halo %d is empty of stellar matter"%i
	else:
	    weights = np.ones(len(xb[i]))
	    normFactor = np.double(len(xb[i].T[0]))
	    x0 = xb[i].T[0] - np.dot(xb[i].T[0],weights)/normFactor
	    x1 = xb[i].T[1] - np.dot(xb[i].T[1],weights)/normFactor
	    x2 = xb[i].T[2] - np.dot(xb[i].T[2],weights)/normFactor
	    tens = np.zeros((3,3))
	    tens[0, 0] = np.dot(x0,x0)/normFactor
	    tens[1, 1] = np.dot(x1,x1)/normFactor        
	    tens[2, 2] = np.dot(x2,x2)/normFactor
	    tens[1, 0] = np.dot(x1,x0)/normFactor
	    tens[2, 0] = np.dot(x2,x0)/normFactor
	    tens[2, 1] = np.dot(x2,x1)/normFactor
	    tens[0, 2] = tens[2, 0]
	    tens[0, 1] = tens[1, 0]
	    tens[1, 2] = tens[2, 1]
	    # Evaluate the mass-weighted centroid along each axis
	    X = np.trapz(mb[i]*xb[i].T[0], xb[i].T[0])/np.trapz(mb[i], xb[i].T[0])
	    Y = np.trapz(mb[i]*xb[i].T[1], xb[i].T[1])/np.trapz(mb[i], xb[i].T[1])
	    Z = np.trapz(mb[i]*xb[i].T[2], xb[i].T[2])/np.trapz(mb[i], xb[i].T[2])
	    # Compute the eigenvalues of the halos
	    w, v = np.linalg.eigh(tens)
	    eigvalues[1,i] = w
	    eigvectors[1,i] = v
	    centroids[1,i] = np.array([X,Y,Z])
    print "Saving output"
    #np.save("eigenvaluesb.npy", eigvalues)
    #np.save("eigenvectorsb.npy", eigvectors)
    #np.save("centroids.npy", centroids)

    out = fi.FITS('subhalo_cat.fits','rw')

    dat=np.zeros(len(h), dtype=[('x', float), ('y', float), ('z', float), ('lambda1_x', float), ('lambda1_y', float), ('lambda1_z', float), ('lambda2_x', float), ('lambda2_y', float), ('lambda2_z', float), ('lambda3_x', float), ('lambda3_y', float), ('lambda3_z', float)])
    dat['lambda1_x'] = eigvectors[0].T[0,0]
    dat['lambda1_y'] = eigvectors[0].T[0,1]
    dat['lambda1_z'] = eigvectors[0].T[0,2]
    dat['lambda2_x'] = eigvectors[0].T[1,0]
    dat['lambda2_y'] = eigvectors[0].T[1,1]
    dat['lambda2_z'] = eigvectors[0].T[1,2]
    dat['lambda3_x'] = eigvectors[0].T[2,0]
    dat['lambda3_y'] = eigvectors[0].T[2,1]
    dat['lambda3_z'] = eigvectors[0].T[2,2]
    dat['x'] = centroids[0].T[0]
    dat['y'] = centroids[0].T[1]
    dat['z'] = centroids[0].T[2]

    out.write(dat)
    out[-1].write_key('EXTNAME', 'baryons')

    dat2=np.zeros(len(h), dtype=[('x', float), ('y', float), ('z', float), ('lambda1_x', float), ('lambda1_y', float), ('lambda1_z', float), ('lambda2_x', float), ('lambda2_y', float), ('lambda2_z', float), ('lambda3_x', float), ('lambda3_y', float), ('lambda3_z', float)])
    dat2['lambda1_x'] = eigvectors[1].T[0,0]      
    dat2['lambda1_y'] = eigvectors[1].T[0,1]
    dat2['lambda1_z'] = eigvectors[1].T[0,2]
    dat2['lambda2_x'] = eigvectors[1].T[1,0]
    dat2['lambda2_y'] = eigvectors[1].T[1,1]
    dat2['lambda2_z'] = eigvectors[1].T[1,2]
    dat2['lambda3_x'] = eigvectors[1].T[2,0]
    dat2['lambda3_y'] = eigvectors[1].T[2,1]
    dat2['lambda3_z'] = eigvectors[1].T[2,2]
    dat2['x'] = centroids[1].T[0]
    dat2['y'] = centroids[1].T[1]
    dat2['z'] = centroids[1].T[2]

    out.write(dat2)
    out[-1].write_key('EXTNAME', 'baryons')
    out.close()

