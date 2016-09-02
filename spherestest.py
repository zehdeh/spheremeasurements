#! bin/python

#from mayavi.mlab import *
import sys
import math, random
from mpl_toolkits.mplot3d import Axes3D, proj3d
import matplotlib.pyplot as plt
from matplotlib.mlab import PCA
from matplotlib.tri import Triangulation
from matplotlib import cm
import numpy as np
import numpy.linalg
from numpy.linalg import norm
from scipy.spatial import ConvexHull
import scipy.interpolate
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle
from mplcolorhelper import MplColorHelper
from opendr.serialization import load_mesh
from opendr.geometry import GaussianCurvature
from utils import *



def plotRadialDeviation(fig,sphere, center, radius, minRadius, maxRadius):
	zValues = sphere[2]-center[2]
	h = np.sqrt((sphere[0]-center[0])**2 + (sphere[2]-center[2])**2)
	r = distance(sphere,center)

	alphas = np.arcsin(zValues/h)
	zIndices = np.where(zValues<0)
	z2Indices = np.where(zValues>0)
	
	xIndices = np.where((sphere[0]-center[0])<0)
	#alphas[zIndices] += np.pi
	alphas[xIndices] = alphas[xIndices]*-1+np.pi

	ax = fig.add_subplot(211, projection='polar')
	ax.scatter(alphas, r)
	ax.grid(True)

def projectToPlane(vertices, pointOnPlane, planeNormalVector):
	vecToCenter = vertices.transpose() - pointOnPlane[0:3]
	distanceToPlane = np.dot(vecToCenter, planeNormalVector)
	diffValues = (planeNormalVector * distanceToPlane[:,np.newaxis]).transpose()[1]
	return (vertices.transpose() - (planeNormalVector * distanceToPlane[:,np.newaxis])).transpose()


def PCA(vertices):
	cov = np.cov(vertices)
	eigval, eigvec = np.linalg.eig(cov)
	return eigval, eigvec.T

def plotEigenVectors(centerPoint, eigval, eigvec, ax):
	for vec,val in zip(eigvec,eigval):
		arrowEndPoint = centerPoint[0:3] + vec*val*0.1
		arrow = Arrow3D([centerPoint[0], arrowEndPoint[0]], [centerPoint[1], arrowEndPoint[1]], [centerPoint[2], arrowEndPoint[2]], mutation_scale=20, lw=2, arrowstyle='-|>', color="r")
		ax.add_artist(arrow)

def projectPoints(vertices, limitsMin, ax):
	for i,axis in zip(range(3),['x','y','z']):
		others = [j for j in range(3) if j != i]
		ax.scatter(vertices[others[0]], vertices[others[1]], zdir=axis, zs=limitsMin[i],cmap=cm.coolwarm)

def plot3D(vertices, triangles, curvature, centerPoint, meanRadius,minRadius, maxRadius):
	maxCurvature = np.mean(curvature[triangles], axis=1)
	maxCurvature = maxCurvature# / np.max(maxCurvature)
	colorHelper = MplColorHelper('coolwarm', np.min(maxCurvature), np.max(maxCurvature), maxCurvature)
	polyColors = [colorHelper.get_rgb(c) for c in maxCurvature]

	poly3dCollection = Poly3DCollection(vertices.T[triangles], facecolors=polyColors, edgecolors='black')


	fig = plt.figure()
	ax3d = fig.add_subplot(221, projection='3d')
	ax3d.add_collection(poly3dCollection)

	means = [np.mean(v) for v in vertices]
	relativeMins = [np.min(v - m) for v,m in zip(vertices,means)]
	relativeMaxs = [np.max(v - m) for v,m in zip(vertices,means)]
	limitsMin = np.add(means,relativeMins)
	limitsMax = np.add(means,relativeMaxs)
	valMin = np.min(np.add(means,relativeMins))
	valMax = np.max(np.add(means,relativeMaxs))
	plt.colorbar(colorHelper.get_mappable())

	axX = fig.add_subplot(222)
	axY = fig.add_subplot(223)
	axZ = fig.add_subplot(224)

	for axisI,axisS,ax in zip(range(3),['x','y','z'], [axX, axY, axZ]):
		otherAxis = [j for j in range(3) if j != axisI]
		circle = Circle((centerPoint[otherAxis[0]], centerPoint[otherAxis[1]]), meanRadius, fill=False, zorder=-1)
		ax.add_patch(circle)

		maxCircle = Circle((centerPoint[otherAxis[0]], centerPoint[otherAxis[1]]), maxRadius, fill=False, color='r', zorder=-1)
		ax.add_patch(maxCircle)

		minCircle = Circle((centerPoint[otherAxis[0]], centerPoint[otherAxis[1]]), minRadius, fill=False, color='b', zorder=-1)
		ax.add_patch(minCircle)

		ax.scatter(vertices[otherAxis[0]], vertices[otherAxis[1]], c=curvature, cmap=cm.coolwarm)


	#projectPoints(vertices, limitsMin, ax3d)

	ax3d.set_xlim(limitsMin[0], limitsMax[0])
	ax3d.set_ylim(limitsMin[1], limitsMax[1])
	ax3d.set_zlim(limitsMin[2], limitsMax[2])

	ax3d.set_xlabel('x')
	ax3d.set_ylabel('y')
	ax3d.set_zlabel('z')
	plt.show()

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def unit_vector(vector):
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def projectByAngleToPlane(vertices, maxAngle, centerPoint, dimensionNo):
	maxAngle = maxAngle * (np.pi/180)
	verticesOnPlane = vertices.copy()
	verticesOnPlane[dimensionNo] = 0


	directionVectors = vertices.T - centerPoint[0:3]
	planeVectors = verticesOnPlane.T - centerPoint[0:3]
	angles = [angle_between(v1, v2) for v1,v2 in zip(directionVectors, planeVectors)]
	print np.degrees(angles)

	indicesBelowMaxAngle = [i for i,a in enumerate(angles) if abs(a) < maxAngle]
	for i in indicesBelowMaxAngle:
		axis = np.cross(directionVectors[i], planeVectors[i])
		if axis[0] == 0 and axis[1] == 0 and axis[2] == 0: continue
		newPoint = np.dot(rotation_matrix(axis,angles[i]),directionVectors[i])
		vertices.T[i] = newPoint + centerPoint[0:3]

if __name__ == '__main__':
	if len(sys.argv) < 2:
		vertices,faces = generateSphere(600, 10, False)
	else:
		mesh = load_mesh(sys.argv[1])
		vertices = mesh.v.T
		faces = mesh.f

	curvature = np.asarray(GaussianCurvature(vertices.T, faces))
	#curvature = curvature / len(vertices.T)
	print 'Curvature min/max: ' + str(np.min(curvature)) + '/' + str(np.max(curvature))
	averageCurvature = np.mean(curvature)
	print 'Average curvature: ' + str(averageCurvature)

	#eigval, eigvec = PCA(vertices)



	meanRadius, centerPoint = fitSphere(vertices)
	print 'Centerpoint: ' + str(centerPoint)
	vertices = vertices.T - centerPoint
	vertices = vertices.T

	centerPoint = [0, 0, 0]

	sphericalCoordinates = toSphericalCoordinates(vertices.T)
	sortedIndices = np.lexsort((sphericalCoordinates.T[2], sphericalCoordinates.T[1]))
	sphericalCoordinates = [[sphericalCoordinates.T[0][i], sphericalCoordinates.T[1][i], sphericalCoordinates.T[2][i]] for i in sortedIndices]
	sphericalMap = np.asarray(sphericalCoordinates)[:,1:3].T
	thetavals,phivals = sphericalMap
	X,Y = np.meshgrid(thetavals,phivals)
	rbf = scipy.interpolate.Rbf(thetavals, phivals, curvature[sortedIndices], function='gaussian')
	Z = rbf(X,Y)

	fig = plt.figure()
	plt.imshow(Z, vmin=Z.min(), vmax=Z.max(), origin='lower', extent=[thetavals.min(), thetavals.max(), phivals.min(), phivals.max()])
	plt.show()

	[minRadiusDev,maxRadiusDev, totalDev] = measureRadialDeviation(vertices, centerPoint, meanRadius)
	expectedCurvature = (1/(meanRadius*meanRadius)) / len(vertices.T)
	print 'Expected curvature: ' + str(expectedCurvature)
	print 'Diff to expected curvature: ' + str(abs(expectedCurvature - averageCurvature))
	print 'Total curvature: ' + str(expectedCurvature*len(vertices.T))

	print 'Expected radius: ' + str(35)
	print 'Fitted radius: ' + str(meanRadius)
	print 'Radial deviation: ' + str(maxRadiusDev+minRadiusDev) + ' (+' + str(maxRadiusDev) + ', -' + str(minRadiusDev) + ')'
	#projectByAngleToPlane(vertices, 45, centerPoint, 2)

	plot3D(vertices, faces, curvature, centerPoint, meanRadius, meanRadius - minRadiusDev, meanRadius + maxRadiusDev)
