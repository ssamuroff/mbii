import numpy as np 
import mbii.lego_tools as utils

def test4():
	"""Reposition the galaxy on the x axis, as defined by the host halo.
	   Set the rotated position to sit on one of the four coordinate axes.
	   The idea is to test the shape rotation function is behaving correctly."""
	phi_shifted = [0, np.pi/2, np.pi]
	theta_shifted = [-np.pi, -np.pi/2,  0, np.pi/2, np.pi]

	pos = np.array([0.5, 0, 0])
	R = np.sqrt(sum(pos*pos))

	data = np.zeros(1, dtype=[('a1',float), ('a2',float), ('a3',float)])
	data['a1'] = 1
	data['a2'] = 0
	data['a3'] = 0

	print 'Starting at '
	print 'position = ', pos
	print 'orientation = ', data

	for p in phi_shifted:
		for t in theta_shifted:
			print '--'
			print 'theta : %3.3f, phi : %3.3f'%(t,p)
			# Galaxy position in radial coordinates
			phi = np.arccos(pos[2]/R) * pos[0]/abs(pos[0])
			theta = np.arcsin(pos[1]/R/np.sin(phi))

			# Work out the rotated galaxy position
			xrot = R * np.sin(p) * np.cos(t)
			yrot = R * np.sin(p) * np.sin(t)
			zrot = R * np.cos(p)

			rotated = np.array([xrot, yrot, zrot])

			# Then do the 3D shape vector
			arot = utils.rotate_shape_vector(data, t , p, theta, phi, axis='a%d')
			print 'position = %3.3f %3.3f %3.3f'%(rotated[0], rotated[1], rotated[2])
			print 'orientation = %3.3f %3.3f %3.3f'%(arot[0], arot[1], arot[2])
			print ''

def test5():
	# Perfect sphere
	a3d = np.array([[1,0,0]])
	b3d = np.array([[0,1,0]])
	c3d = np.array([[0,0,1]])

	q3d = np.array([1])
	s3d = np.array([1])

	e1,e2 = utils.project_3d_shape(a3d, b3d, c3d, q3d, s3d)
	print 'Spherical'
	print 'e1,e2 = ', e1[0], e2[0]
	print '--'


	a3d = np.array([[0,0,1]])
	b3d = np.array([[1,0,0]])
	c3d = np.array([[0,1,0]])

	q3d = np.array([0.5])
	s3d = np.array([0.5])

	e1,e2 = utils.project_3d_shape(a3d, b3d, c3d, q3d, s3d)
	print 'Oblate along the z axis'
	print 'e1,e2 = ', e1[0], e2[0]


	a3d = np.array([[1,0,0]])
	b3d = np.array([[0,0,1]])
	c3d = np.array([[0,1,0]])

	q3d = np.array([0.5])
	s3d = np.array([0.5])

	e1,e2 = utils.project_3d_shape(a3d, b3d, c3d, q3d, s3d)
	print 'Oblate along the x axis'
	print 'e1,e2 = ', e1[0], e2[0]

	a3d = np.array([[0,1,0]])
	b3d = np.array([[0,0,1]])
	c3d = np.array([[1,0,0]])

	q3d = np.array([0.5])
	s3d = np.array([0.5])

	e1,e2 = utils.project_3d_shape(a3d, b3d, c3d, q3d, s3d)
	print 'Oblate along the y axis'
	print 'e1,e2 = ', e1[0], e2[0]
	print '--'

	a3d = np.array([[0,0,1]])
	b3d = np.array([[0.5,0.5,0]])
	c3d = np.array([[0.5,-0.5,0]])

	q3d = np.array([0.5])
	s3d = np.array([0.25])

	e1,e2 = utils.project_3d_shape(a3d, b3d, c3d, q3d, s3d)
	print 'Triaxial, longer dimension along the z axis, mid dimension at 45 degrees to the x axis'
	print 'e1,e2 = ', e1[0], e2[0]






