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

def MTF(x):
	#y = np.diff(x)

	#y = np.append(x, np.zeros(100))
	Y = np.fft.fft(x)

	return Y#[:len(Y) // 2]

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

def getPlaneIntersection(n1, n2, c1, c2):
	u = np.cross(n1,n2)
	maxc = np.argmax(np.fabs(n1))

	d1 = -n1.dot(c1)
	d2 = -n2.dot(c2)

	iP = np.zeros(3)

	if maxc == 0:
		iP[1] = (d2*n1[2] - d1*n2[2]) / u[0]
		iP[2] = (d1*n2[1] - d2*n1[1]) / u[0]
	elif maxc == 1:
		iP[0] = (d1*n2[2] - d2*n1[2]) / u[1]
		iP[2] = (d2*n1[0] - d1*n2[0]) / u[1]
	else:
		iP[0] = (d2*n1[1] - d1*n2[1]) / u[2]
		iP[1] = (d1*n2[0] - d2*n1[0]) / u[2]
	
	return iP, iP+1000*u

def plotProjection(v, vProjected, n1, n2, plane1, plane2, c1, c2, centerPoint, crossProduct, upVector, xVector, ip, perp1, perp2):
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	[xx, yy] = np.meshgrid(range(-80,20),range(-50,50))
	z = (-plane1[0]*xx - plane1[1]*yy - plane1[3])/plane1[2]
	#ax.plot_surface(xx,yy,z)

	p0 = c1
	arrow = Arrow3D(\
	[p0[0],p0[0]+n1[0]*50],\
	[p0[1],p0[1]+n1[1]*50],\
	[p0[2],p0[2]+n1[2]*50],\
	mutation_scale=20,\
	arrowstyle="-|>")
	ax.add_artist(arrow)

	z = (-plane2[0]*xx - plane2[1]*yy - plane2[3])/plane2[2]
	#ax.plot_surface(xx,yy,z,color='r')

	p0 = c2
	arrow = Arrow3D(\
	[p0[0],p0[0]+n2[0]*50],\
	[p0[1],p0[1]+n2[1]*50],\
	[p0[2],p0[2]+n2[2]*50],\
	mutation_scale=20,\
	arrowstyle="-|>")
	ax.add_artist(arrow)

	p0 = ip
	arrow = Arrow3D(\
	[p0[0],p0[0]+crossProduct[0]*50],\
	[p0[1],p0[1]+crossProduct[1]*50],\
	[p0[2],p0[2]+crossProduct[2]*50],\
	mutation_scale=20,\
	arrowstyle="-|>")
	ax.add_artist(arrow)

	arrow = Arrow3D(\
	[p0[0],p0[0]+upVector[0]*50],\
	[p0[1],p0[1]+upVector[1]*50],\
	[p0[2],p0[2]+upVector[2]*50],\
	mutation_scale=20,\
	arrowstyle="-|>")
	ax.add_artist(arrow)

	arrow = Arrow3D(\
	[p0[0],p0[0]+xVector[0]*50],\
	[p0[1],p0[1]+xVector[1]*50],\
	[p0[2],p0[2]+xVector[2]*50],\
	mutation_scale=20,\
	arrowstyle="-|>")
	ax.add_artist(arrow)

	arrow = Arrow3D(\
	[p0[0],p0[0]+perp1[0]*50],\
	[p0[1],p0[1]+perp1[1]*50],\
	[p0[2],p0[2]+perp1[2]*50],\
	color='r',\
	mutation_scale=20,\
	arrowstyle="-|>")
	ax.add_artist(arrow)

	arrow = Arrow3D(\
	[p0[0],p0[0]+perp2[0]*50],\
	[p0[1],p0[1]+perp2[1]*50],\
	[p0[2],p0[2]+perp2[2]*50],\
	color='r',\
	mutation_scale=20,\
	arrowstyle="-|>")
	ax.add_artist(arrow)

	ax.scatter(v.T[0],v.T[1],v.T[2], color='b', marker='.')
	#ax.scatter(vProjected.T[0],vProjected.T[1],vProjected.T[2], color='yellow', marker='o')
	ax.scatter(ip[0],ip[1],ip[2], color='black', marker='o')
	plt.xlabel('x')
	plt.ylabel('y')
	
	plt.show()

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

	verticesLeftCenter = np.mean(verticesLeft, axis=0)
	verticesRightCenter = np.mean(verticesRight, axis=0)
	
	plane1 = fitPlane(verticesLeft.T, verticesLeftCenter.tolist() + [1])
	plane2 = fitPlane(verticesRight.T, verticesRightCenter.tolist() + [1])

	normal1 = plane1[0:3]/np.linalg.norm(plane1[0:3])
	normal2 = plane2[0:3]/np.linalg.norm(plane2[0:3])


	intersectionPoint1,ip2 = getPlaneIntersection(normal1, normal2, verticesLeftCenter[0:3], verticesRightCenter[0:3])

	crossProduct = np.cross(normal1,normal2)
	crossProduct = crossProduct/np.linalg.norm(crossProduct)

	upVector = normal1 + normal2
	upVector = upVector/np.linalg.norm(upVector)

	xVector = np.cross(upVector, crossProduct)
	xVector = xVector/np.linalg.norm(xVector)

	distances = distance.pdist(vertices)
	distances = distance.squareform(distances)
	np.fill_diagonal(distances, np.inf)
	nyquistRate = distances.min(axis=0).min()/4
	nyquistFrequency = 1/(2*nyquistRate)

	verticesUnprojected = vertices
	v = vertices - centerPoint
	dist = v.dot(crossProduct)
	vertices = vertices - dist[...,None]*crossProduct

	intersectionPoint1 = intersectionPoint1 - ((intersectionPoint1 - centerPoint).dot(crossProduct))*crossProduct

	perp1 = np.cross(normal1, crossProduct)
	perp2 = -np.cross(normal2, crossProduct)

	
	plotProjection(verticesUnprojected, vertices, normal1, normal2, plane1, plane2,\
	verticesLeftCenter, verticesRightCenter, centerPoint,\
	crossProduct, upVector, xVector, intersectionPoint1, perp1, perp2)


	xCoordinates = xVector.dot(vertices.T - centerPoint[...,None])
	yCoordinates = upVector.dot(vertices.T - centerPoint[...,None])
	intersectionPoint = [xVector.dot(intersectionPoint1 - centerPoint), upVector.dot(intersectionPoint1 - centerPoint)]
	perp1 = [xVector.dot(perp1), upVector.dot(perp1)]
	perp2 = [xVector.dot(perp2), upVector.dot(perp2)]

	perp1 = perp1/np.linalg.norm(perp1)
	perp2 = perp2/np.linalg.norm(perp2)

	angularFactor1 = 1/np.array([1,0]).dot(perp1)*np.linalg.norm(perp1)
	angularFactor2 = 1/np.array([-1,0]).dot(perp2)*np.linalg.norm(perp2)

	xShift = math.fabs(xCoordinates.min())
	xCoordinates += xShift
	intersectionPoint[0] += xShift

	xSpread = xCoordinates.max()
	binSize = nyquistRate
	noBins = int(pow(2, math.floor(math.log(xSpread / binSize) / math.log(2))))

	perfectProfile1 = np.array([intersectionPoint + (binSize/2)*angularFactor1*perp1 + i*binSize*angularFactor1*perp1 for i in range(0,noBins/2)]).T
	perfectProfile2 = np.array([intersectionPoint - (binSize/2)*angularFactor2*perp2 + i*binSize*angularFactor2*perp2 for i in range((noBins/2),0,-1)]).T

	binsPerfectProfile = np.hstack((perfectProfile2[1], perfectProfile1[1]))
	xTruePeak = intersectionPoint[0]

	binsMeasuredProfile = np.zeros(noBins)

	fig = plt.figure(1)
	plt.subplot(211)
	plt.plot([xTruePeak, xTruePeak], [yCoordinates.min(),intersectionPoint[1]], color='black')
	plt.scatter(xCoordinates, yCoordinates)

	for i in range(0, noBins/2):
		xMinLeft = xTruePeak - (i+1)*binSize
		xMaxLeft = xTruePeak - i*binSize
		plt.plot([xMinLeft, xMinLeft], [yCoordinates.min(),intersectionPoint[1]], color='lightgray')
		yInBinLeft = yCoordinates[np.all([xMinLeft < xCoordinates,xCoordinates < xMaxLeft], axis=0)]
		if len(yInBinLeft) > 0:
			binsMeasuredProfile[noBins/2 - i - 1] = yInBinLeft.mean()
			plt.scatter(xMinLeft + (xMaxLeft - xMinLeft)/2,binsMeasuredProfile[noBins/2 - i - 1], color='g')

		xMinRight = xTruePeak + i*binSize
		xMaxRight = xTruePeak + (i+1)*binSize
		plt.plot([xMaxRight, xMaxRight], [yCoordinates.min(),intersectionPoint[1]], color='lightgray')
		yInBinRight = yCoordinates[np.all([xMinRight < xCoordinates,xCoordinates < xMaxRight], axis=0)]
		if len(yInBinRight) > 0:
			binsMeasuredProfile[noBins / 2 + i] = yInBinRight.mean()
			plt.scatter(xMinRight + (xMaxRight - xMinRight)/2,binsMeasuredProfile[noBins/2 + i], color='g')

	plt.plot([xTruePeak-(noBins/2)*binSize, xTruePeak-(noBins/2)*binSize], [yCoordinates.min(),intersectionPoint[1]], color='black')
	plt.plot([xTruePeak+(noBins/2)*binSize, xTruePeak+(noBins/2)*binSize], [yCoordinates.min(),intersectionPoint[1]], color='black')
	
	for i in range(noBins):
		if binsMeasuredProfile[i] == 0 and i+1 < noBins:
			binsMeasuredProfile[i] = (binsMeasuredProfile[i-1] + binsMeasuredProfile[i+1])/2
	
	plt.scatter(intersectionPoint[0], intersectionPoint[1], color='black')
	#plt.scatter(perfectProfile1[0], perfectProfile1[1], color='r')
	#plt.scatter(perfectProfile2[0], perfectProfile2[1], color='r')

	plt.subplot(212)

	plt.scatter(np.arange(0, noBins), binsMeasuredProfile)
	plt.scatter(np.arange(0, noBins), binsPerfectProfile, color='r')

	xShift = min(binsMeasuredProfile)
	binsMeasuredProfile = binsMeasuredProfile - xShift
	binsPerfectProfile = binsPerfectProfile - xShift

	binsMeasuredProfile[0] = 0
	binsMeasuredProfile[-1] = 0
	binsPerfectProfile[0] = 0
	binsPerfectProfile[-1] = 0

	j = np.arange(0,noBins)
	wj = 1 - ((j - noBins/2.)/(noBins/2.))**2
	binsMeasuredProfile = binsMeasuredProfile*wj
	binsPerfectProfile = binsPerfectProfile*wj
	
	binsMeasuredProfile = np.append(binsMeasuredProfile, np.zeros(noBins))
	binsPerfectProfile = np.append(binsPerfectProfile, np.zeros(noBins))
	for i in range(noBins):
		binsMeasuredProfile[noBins*2 - i - 1] = -binsMeasuredProfile[i]
		binsPerfectProfile[noBins*2 - i - 1] = -binsPerfectProfile[i]

	plt.show()
	fig = plt.figure()
	plt.subplot(211)

	plt.plot(np.arange(0,binsMeasuredProfile.shape[0]), binsMeasuredProfile)
	plt.plot(np.arange(0,binsPerfectProfile.shape[0]), binsPerfectProfile, linestyle='--')

	plt.subplot(212)

	Ymeasured = np.fft.fft(binsMeasuredProfile)
	Yperfect = np.fft.fft(binsPerfectProfile)

	H = (Ymeasured.imag[1::2]/Yperfect.imag[1::2])
	frequencies = np.fft.fftfreq(H.shape[0],d=nyquistRate/2)

	freqIdxAboveZero = frequencies >= 0
	frequencies = frequencies[freqIdxAboveZero]
	H = H[freqIdxAboveZero]
	#frequencies = np.fft.fftshift(frequencies)

	plt.xlim([0,nyquistFrequency*1.2])
	plt.ylim([0,1.1])
	plt.plot([nyquistFrequency,nyquistFrequency], [H.min(),H.max()])
	plt.plot(frequencies, H)

	#plt.plot(np.arange(0,Ymeasured.shape[0]), Ymeasured)
	#plt.plot(np.arange(0,Yperfect.shape[0]), Yperfect, linestyle='--')


	#plt.plot(np.arange(0, Ymeasured[::2].shape[0]), np.fabs(Ymeasured[::2].real))
	#plt.plot(np.arange(0, Yperfect[::2].shape[0]), np.fabs(Yperfect[::2].real), linestyle='--')
	#plt.plot(np.fft.fftfreq(Ymeasured.shape[0], d=nyquistRate/2), abs(Ymeasured))
	#plt.plot(np.fft.fftfreq(Yperfect.shape[0], d=nyquistRate/2), abs(Yperfect), linestyle='--')

	plt.show()

