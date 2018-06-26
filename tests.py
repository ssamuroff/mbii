import numpy as np 
import mbii.lego_tools as utils
import mbii.symmetrise_lib as lib

def test_rotation_matrix():
	# Start on the x axis, at unit distance
	position = np.array([1,0,0])
	answers = [[1,0,0], [0,1,0], [-1,0,0], [0,-1,0], [-1,0,0] ]
	rotations = [0, np.pi/2, np.pi, -np.pi/2, -np.pi]
	print 'Initial position: ', position
	print ''
	print 'Rotation axis: (theta,phi) = 0,0'

	for rot,expected in zip(rotations, answers):
		Rxyz = lib.build_rotation_matrix(alpha=rot, theta=0.0, phi=0.0)
		rotated = np.dot(Rxyz,position)
		print 'rotation angle: %3.3f rad'%rot,
		if np.isclose(rotated,expected).all():
			print 'Okay.'
		else:
			print 'Not okay.'
			print 'Rotated position does not match the expected one'
			print rotated, expected

	# Try the same thing, but with a rotation vector aligned with the x axis
	answers = [[1,0,0], [1,0,0], [1,0,0], [1,0,0], [1,0,0] ]
	rotations = [0, np.pi/2, np.pi, -np.pi/2, -np.pi]
	print ''
	print 'Rotation axis: (theta,phi) = 0,pi/2'

	for rot,expected in zip(rotations, answers):
		Rxyz = lib.build_rotation_matrix(alpha=rot, theta=np.pi/2, phi=0.0)
		rotated = np.dot(Rxyz,position)
		print 'rotation angle: %3.3f rad'%rot,
		if np.isclose(rotated,expected).all():
			print 'Okay.'
		else:
			print 'Not okay.'
			print 'Rotated position does not match the expected one'
			print rotated, expected

	# Try the same thing, but with a rotation vector aligned with the (negative) y axis
	answers = [[1,0,0], [0,-1,0], [-1,0,0], [0,1,0], [-1,0,0] ]
	rotations = [0, np.pi/2, np.pi, -np.pi/2, -np.pi]
	print ''
	print 'Rotation axis: (theta,phi) = 0,pi'

	for rot,expected in zip(rotations, answers):
		Rxyz = lib.build_rotation_matrix(alpha=rot, theta=np.pi, phi=0.0)
		rotated = np.dot(Rxyz,position)
		print 'rotation angle: %3.3f rad'%rot,
		if np.isclose(rotated,expected).all():
			print 'Okay.'
		else:
			print 'Not okay.'
			print 'Rotated position does not match the expected one'
			print rotated, expected

def test_sphericity(cat, cat0=[]):
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.pyplot as plt

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	mask = cat['most_massive']==0
	ax.plot(cat['x'][mask], cat['y'][mask], cat['z'][mask], '.', color='purple', alpha=0.5)
	ax.plot(cat['x'][np.invert(mask)], cat['y'][np.invert(mask)], cat['z'][np.invert(mask)], 'x', color='pink')

	if len(cat0)>0:
		ax.plot(cat0['x'][mask], cat0['y'][mask], cat0['z'][mask], '*', color='hotpink')

	plt.savefig('/home/ssamurof/spherical_test.png')

	# Work out the quadrupole moments of the distribution
	x = cat['x'][mask]
	y = cat['y'][mask]
	z = cat['z'][mask]
	x0 = cat['x'][mask].mean()
	y0 = cat['y'][mask].mean()
	z0 = cat['z'][mask].mean()
	r2 = (x-x0)**2+(y-y0)**2+(z-z0)**2
	Qxx = sum(3*(x-x0)*(x-x0) - r2)
	Qyy = sum(3*(y-y0)*(y-y0) - r2)
	Qzz = sum(3*(z-z0)*(z-z0) - r2)

	# Do the same with a trefoil configuration
	nelem = len(x)
	r0 = (x-x0).max()
	pad = np.zeros(nelem)
	dx=[]
	dy=[]
	dz=[]
	for i in xrange(nelem):
		r = np.random.choice([r0,-r0])
		j = np.random.choice(3)
		if j==0:
			dx.append(r) ; dy.append(0) ; dz.append(0)
		if j==1:
			dx.append(0) ; dy.append(r) ; dz.append(0)
		if j==2:
			dx.append(0) ; dy.append(0) ; dz.append(r)

	dx = np.array(dx)
	dy = np.array(dy)
	dz = np.array(dz)

	r2 = dx**2+dy**2+dz**2
	Qxx_ref = sum(3*dx*dx - r2)
	Qyy_ref = sum(3*dy*dy - r2)
	Qzz_ref = sum(3*dz*dz - r2)

	print Qxx/Qzz,Qyy/Qzz,Qzz/Qzz
	print 'Reference case: ', Qxx_ref/Qzz_ref,Qyy_ref/Qzz_ref,Qzz_ref/Qzz_ref







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

def test_projection():
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






