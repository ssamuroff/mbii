import numpy as np
from readsubhalo import *

def compute_inertia_tensors(snap):
    """ Compute the intertia tensors for all subhalos for
    the dark matter and stellar components and saves the output for
    latter integration within the database
    """
    h = snap.readsubhalo()
    x = snap.load(1, 'pos', h)
    xb = snap.load(4, 'pos', h)
    
    eigvectors = np.zeros((2, len(h), 3, 3))
    eigvalues  = np.zeros((2, len(h), 3))
 
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
	    # Compute the eigenvalues of the halos
	    w, v = np.linalg.eigh(tens)
	    eigvalues[0,i] = w
	    eigvectors[0,i] = v
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
	    # Compute the eigenvalues of the halos
	    w, v = np.linalg.eigh(tens)
	    eigvalues[1,i] = w
	    eigvectors[1,i] = v
    print "Saving output"
    np.save("eigenvaluesb.npy", eigvalues)
    np.save("eigenvectorsb.npy", eigvectors)
    
