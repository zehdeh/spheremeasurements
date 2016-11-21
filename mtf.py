#! /usr/bin/python2.7

import sys
import math
import numpy as np
from src.OBJIO import loadOBJ
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from src.utils import Arrow3D
from src.fitting import fitPlane
from scipy.spatial import distance

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
	
	vertices, faces, normals = loadOBJ(sys.argv[1])
	centerPoint = np.mean(vertices, axis=0)
	vertices = vertices - centerPoint

	indicesLeft = np.where(vertices.T[2] > 15)
	indicesRight = np.where(vertices.T[2] < -5)

	verticesLeft = vertices[indicesLeft]
	verticesRight = vertices[indicesRight]

	verticesLeftCenter = np.mean(verticesLeft, axis=0).tolist() + [1]
	verticesRightCenter = np.mean(verticesRight, axis=0).tolist() + [1]
	
	plane1 = fitPlane(verticesLeft.T, verticesLeftCenter)
	plane2 = fitPlane(verticesRight.T, verticesRightCenter)

	fig = plt.figure()
	#ax = fig.add_subplot(111, projection='3d')

	[xx, yy] = np.meshgrid(range(100),range(100))
	z = (-plane1[0]*xx - plane1[1]*yy - plane1[3])/plane1[2]
	#ax.plot_surface(xx,yy,z)

	p0 = verticesLeftCenter
	normal1 = plane1[0:3]/np.linalg.norm(plane1[0:3])
	arrow = Arrow3D(\
	[p0[0],p0[0]+normal1[0]*100],\
	[p0[1],p0[1]+normal1[1]*100],\
	[p0[2],p0[2]+normal1[2]*100])
	#ax.add_artist(arrow)

	z = (-plane2[0]*xx - plane2[1]*yy - plane2[3])/plane2[2]
	#ax.plot_surface(xx,yy,z,color='r')

	p0 = verticesRightCenter
	normal2 = plane2[0:3]/np.linalg.norm(plane2[0:3])
	arrow = Arrow3D(\
	[p0[0],p0[0]+normal2[0]*100],\
	[p0[1],p0[1]+normal2[1]*100],\
	[p0[2],p0[2]+normal2[2]*100])
	#ax.add_artist(arrow)

	crossProduct = np.cross(normal1,normal2)
	crossProduct = crossProduct/np.linalg.norm(crossProduct)

	p0 = centerPoint
	arrow = Arrow3D(\
	[p0[0],p0[0]+crossProduct[0]*100],\
	[p0[1],p0[1]+crossProduct[1]*100],\
	[p0[2],p0[2]+crossProduct[2]*100])
	#ax.add_artist(arrow)

	upVector = normal1 + normal2
	upVector = upVector/np.linalg.norm(upVector)
	arrow = Arrow3D(\
	[p0[0],p0[0]+upVector[0]*100],\
	[p0[1],p0[1]+upVector[1]*100],\
	[p0[2],p0[2]+upVector[2]*100])
	#ax.add_artist(arrow)

	xVector = np.cross(upVector, crossProduct)
	xVector = xVector/np.linalg.norm(xVector)
	arrow = Arrow3D(\
	[p0[0],p0[0]+xVector[0]*100],\
	[p0[1],p0[1]+xVector[1]*100],\
	[p0[2],p0[2]+xVector[2]*100])
	#ax.add_artist(arrow)

	v = vertices - centerPoint
	dist = v.dot(crossProduct)
	vertices = vertices - dist[...,None]*crossProduct


	xCoordinates = xVector.dot(vertices.T - centerPoint[...,None])
	yCoordinates = upVector.dot(vertices.T - centerPoint[...,None])

	xCoordinates += math.fabs(xCoordinates.min())


	distances = distance.pdist(vertices)
	distances = distance.squareform(distances)
	np.fill_diagonal(distances, np.inf)
	minimalDistance = distances.min()

	xSortedIndices = np.argsort(xCoordinates)
	yPeak = np.argsort(yCoordinates)[-1]
	leftLow = yCoordinates[xSortedIndices][:yPeak].min()
	rightLow = yCoordinates[xSortedIndices][yPeak:].min()
	yCutOff = max(leftLow, rightLow)
	yCoordinates += math.fabs(yCutOff)


	includedIndices = np.where(yCoordinates >= 0)

	xCoordinates = xCoordinates[includedIndices]
	yCoordinates = yCoordinates[includedIndices]




	xSpread = xCoordinates.max()
	noBinsUnadjusted = xSpread / (minimalDistance/2)

	N = xCoordinates.shape[0]
	noBins = int(pow(2, math.ceil(math.log(N)/math.log(2)-5)))

	xTruePeak = xCoordinates[np.argsort(yCoordinates)[-10]].mean()
	print xTruePeak

	bins = np.zeros(noBins)
	binSize = xSpread / noBins
	for i in range(noBins):
		xMin = max(0, (i-1)*binSize)
		xMax = i*binSize
		yInBin = yCoordinates[np.all([xMin < xCoordinates,xCoordinates < xMax], axis=0)]
		if yInBin.shape[0] > 0:
			bins[i] = yInBin.mean()

	plt.plot(np.arange(0, bins.shape[0]), bins)
	#plt.scatter(xCoordinates, yCoordinates)


	#ax.scatter(vertices.T[0],vertices.T[1],vertices.T[2], color='b', marker='.')

	#ax.scatter(verticesLeft.T[0],verticesLeft.T[1],verticesLeft.T[2], color='b', marker='.')
	#ax.scatter(verticesRight.T[0],verticesRight.T[1],verticesRight.T[2], color='r', marker='.')

	plt.show()

'''
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
	#vertices -= distances*normal

	rotationAxis = np.cross([0,0,1], normal)
	rotationAngle = np.arccos(np.dot([0,0,1], normal))
	R = getRotationMatrix(rotationAngle, rotationAxis)
	app = np.ones(vertices.shape[0])
	#vertices = np.vstack([vertices.T, app]).T
	#vertices = R.dot(vertices.T).T
	#vertices[:,2] = 0

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
	#vertices = R.dot(vertices.T)[0:2].T
	#vertices = vertices[vertices[:,0].argsort()].T
	#maxY = max(vertices[1,-1], vertices[1,0])
	#vertices = vertices.T[vertices[1] >= maxY]



	#arrow2 = Arrow3D(\
	#[p0[0], p0[0]+n2[0]*100],\
	#[p0[1], p0[1]+n2[1]*100],\
	#[p0[2], p0[2]+n2[2]*100])

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	#ax = fig.add_subplot(111)
	#ax.plot_surface(xx, yy, z1)
	#ax.plot_surface(xx, yy, z2)
	#ax.add_artist(arrow)
	#ax.add_artist(arrow2)

	ax.scatter(vertices.T[0],vertices.T[1],vertices.T[2], color='b', marker='.')
	#ax.scatter(vertices.T[0],vertices.T[1], color='b', marker='.')
	#ax.plot(vertices.T[0], vertices.T[1])

	#ax.scatter(centerPoint1.T[0],centerPoint1.T[1],centerPoint1.T[2], color='g', marker='o')
	#ax.scatter(centerPoint2.T[0],centerPoint2.T[1],centerPoint2.T[2], color='g', marker='o')
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	#ax.set_zlabel('Z Label')
'''

'''
	# Create cubic bounding box to simulate equal aspect ratio
	max_range = np.array([vertices.T[0].max()-vertices.T[0].min(), vertices.T[1].max()-vertices.T[1].min(), vertices.T[2].max()-vertices.T[2].min()]).max()
	Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(vertices.T[0].max()+vertices.T[0].min())
	Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(vertices.T[1].max()+vertices.T[1].min())
	Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(vertices.T[2].max()+vertices.T[2].min())
	# Comment or uncomment following both lines to test the fake bounding box:
	for xb, yb, zb in zip(Xb, Yb, Zb):
		ax.plot([xb], [yb], [zb], 'w')

	
	plt.show()
'''
