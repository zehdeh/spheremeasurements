#! /usr/bin/python2.7

import sys
import math
import numpy as np
from src.OBJIO import loadOBJviaVTK
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from src.utils import Arrow3D
from src.fitting import fitPlane

def getRotationMatrix(angle, axis):
	sina = math.sin(angle)
	cosa = math.cos(angle)
	R = np.diag([cosa, cosa, cosa])
	R += np.outer(axis, axis) * (1.0 - cosa)
	axis *= sina
	R += np.array([[0.0, -axis[2], axis[1]],
	[axis[2], 0.0, -axis[0]],
	[-axis[1], axis[0], 0.0]])
	M = np.identity(4)
	M[:3, :3] = R
	return M

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print 'Too few arguments'
		sys.exit(0)
	
	vertices, faces, normals, polyData = loadOBJviaVTK(sys.argv[1])
	centerPoint = np.mean(vertices, axis=0)
	vertices = vertices - centerPoint

	centerPoint = np.mean(vertices, axis=0)

	noRows = vertices.shape[0]
	p = (np.ones((noRows,1)))
	AB = np.hstack([vertices - centerPoint, p])
	uu, dd, vv = np.linalg.svd(AB)
	B = vv[1,:]/np.linalg.norm(vv[1,0:3])

	normal = B[0:3]

	p0 = centerPoint

	arrow1 = Arrow3D(\
	[p0[0], p0[0]+normal[0]*100],\
	[p0[1], p0[1]+normal[1]*100],\
	[p0[2], p0[2]+normal[2]*100])

	verticesDirections = vertices - centerPoint
	distances = verticesDirections.dot(normal)[...,None]
	vertices -= distances*normal

	rotationAxis = np.cross([0,0,1], normal)
	rotationAngle = np.arccos(np.dot([0,0,1], normal))
	R = getRotationMatrix(rotationAngle, rotationAxis)
	app = np.ones(vertices.shape[0])
	vertices = np.vstack([vertices.T, app]).T
	vertices = R.dot(vertices.T).T
	vertices[:,2] = 0

	centerPoint = np.mean(vertices, axis=0)

	noRows = vertices.shape[0]
	p = (np.ones((noRows,1)))
	AB = np.hstack([vertices - centerPoint])
	uu, dd, vv = np.linalg.svd(AB)
	B = vv[1,:]/np.linalg.norm(vv[1,0:3])

	normal = B[0:3]
	p0 = centerPoint

	rotationAxis = np.cross([0,-1,0], normal)
	rotationAngle = np.arccos(np.dot([0,-1,0], normal))
	R = getRotationMatrix(rotationAngle, rotationAxis)
	vertices = R.dot(vertices.T)[0:2].T
	vertices = vertices[vertices[:,0].argsort()].T
	maxY = max(vertices[1,-1], vertices[1,0])
	vertices = vertices.T[vertices[1] >= maxY]



	#arrow2 = Arrow3D(\
	#[p0[0], p0[0]+n2[0]*100],\
	#[p0[1], p0[1]+n2[1]*100],\
	#[p0[2], p0[2]+n2[2]*100])

	fig = plt.figure()
	#ax = fig.add_subplot(111, projection='3d')
	ax = fig.add_subplot(111)
	#ax.plot_surface(xx, yy, z1)
	#ax.plot_surface(xx, yy, z2)
	#ax.add_artist(arrow)
	#ax.add_artist(arrow2)

	#ax.scatter(vertices.T[0],vertices.T[1],vertices.T[2], color='b', marker='.')
	#ax.scatter(vertices.T[0],vertices.T[1], color='b', marker='.')
	ax.plot(vertices.T[0], vertices.T[1])

	#ax.scatter(centerPoint1.T[0],centerPoint1.T[1],centerPoint1.T[2], color='g', marker='o')
	#ax.scatter(centerPoint2.T[0],centerPoint2.T[1],centerPoint2.T[2], color='g', marker='o')
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	#ax.set_zlabel('Z Label')

	'''
	# Create cubic bounding box to simulate equal aspect ratio
	max_range = np.array([vertices.T[0].max()-vertices.T[0].min(), vertices.T[1].max()-vertices.T[1].min(), vertices.T[2].max()-vertices.T[2].min()]).max()
	Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(vertices.T[0].max()+vertices.T[0].min())
	Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(vertices.T[1].max()+vertices.T[1].min())
	Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(vertices.T[2].max()+vertices.T[2].min())
	# Comment or uncomment following both lines to test the fake bounding box:
	for xb, yb, zb in zip(Xb, Yb, Zb):
		ax.plot([xb], [yb], [zb], 'w')
	'''

	
	plt.show()
