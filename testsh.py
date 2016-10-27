#! /usr/bin/python
USE_CACHE = False

import os
import sys
import timeit
import numpy as np
import time

from math import isnan, sin, degrees
from scipy.special import sph_harm
from src.OBJIO import loadOBJ
from src.fitting import fitSphere
from src.mesh import Mesh
from src.voronoi import getVoronoiArea
from scipy.optimize import leastsq

from cycler import cycler
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as mplcm

from pyshtools.spectralanalysis import SHPowerSpectrum, SHPowerSpectrumDensity
from pyshtools.expand import SHExpandLSQ

def getSphericalCoordinates(vertices, centerPoint):
	vertices = vertices - centerPoint
	r = np.sqrt(vertices.T[0]**2 + vertices.T[1]**2 + vertices.T[2]**2)
	theta = np.arccos(vertices.T[2] / r)
	phi = np.arctan(vertices.T[1]/vertices.T[0])
	return np.array([phi, theta, r])

def getCartesianCoordinates(theta, phi, r, centerPoint):
	x = r*np.sin(theta)*np.cos(phi)
	y = r*np.sin(theta)*np.sin(phi)
	z = r*np.cos(theta)

	coords = np.array([x,y,z])
	coords = coords.T + centerPoint

	return coords

def residuals(aa, Yval, rs):
	return (Yval.dot(aa) - rs)

def complexResiduals(aa, Yval, rs):
	diff = (Yval.dot(aa) - rs)
	return diff.real**2 + diff.imag**2

def getAMatrix(theta, phi, l, filePath, oldA=None):
	n = len(phi)
	
	path = os.path.split(filePath)
	fileName = path[-1][:-4]

	matrixName = fileName + '_' + str(l)
	matrixPath = os.path.join(path[0],matrixName)
	if os.path.isfile(matrixPath + '.npy') and USE_CACHE:
		A = np.load(matrixPath + '.npy')
		print 'Using cache for ' + str(l)
	else:
		A = np.zeros((n, 2*l+1))
		if oldA is None:
			for i in range(n):
				A[i] = [sph_harm(j-l, l, theta[i], phi[i]).real for j in range(0, 2*l+1)]
			#A = np.nan_to_num(A)
		else:
			A[:,1:2*l] = oldA
			for i in range(n):
				A[i,0] = sph_harm(-l, l, theta[i], phi[i]).real
				A[i,-1] = sph_harm(l, l, theta[i], phi[i]).real
		np.save(matrixPath,A)
	
	return A
def sirf(sphericalCoordinates, Lmax, fileName):
	phi, theta, r = sphericalCoordinates

	s = 10
	# areas auf beide Seiten!
	A = getAMatrix(theta, phi, 0, fileName)
	for i in range(1,s):
		A = np.concatenate((A,getAMatrix(theta, phi, i, fileName)), axis=1)
	
	b = np.zeros((Lmax+1)**2)
	b[0:A.shape[1]], flag = leastsq(residuals, b[0:A.shape[1]], args=(A,r))

	np.set_printoptions(threshold=np.inf)
	res = r - A.dot(b[0:A.shape[1]])
	Al = getAMatrix(theta, phi, s-1, fileName)
	for l in range(s, Lmax+1):
		Al = getAMatrix(theta, phi, l, fileName, Al)
		lowerBound = l**2
		upperBound = (l+1)**2

		b[lowerBound:upperBound], flag = leastsq(residuals, b[lowerBound:upperBound], args=(Al,res))
		res = res - Al.dot(b[lowerBound:upperBound])
		A = np.concatenate((A,Al), axis=1)

	ys = np.zeros(Lmax+1)
	for i in range(Lmax+1):
		lowerBound = i**2
		upperBound = (i+1)**2
		ys[i] = np.linalg.norm(b[lowerBound:upperBound], ord=2)
	
	return b, ys, A

def nlsf(sphericalCoordinates, LMax):
	phi, theta, r = sphericalCoordinates

	n = len(phi)
	k = (Lmax + 1)**2
	Y = np.zeros((n,k))
	
	kToLM = np.zeros((k,2))
	for i in range(n):
		for l in range(Lmax+1):
			for m in range(-l,l+1):
				j = l**2 + l + m + 1
				kToLM[j-1,:] = [l,m]
				s = sph_harm(l, m, theta[i], phi[i]).real
				#print 'l ' + str(l)
				#print 'm ' + str(m)
				if not isnan(s):
					Y[i,j-1] = s

	a0 = np.zeros(k)
	a1, flag = leastsq(residuals, a0, args=(Y,r))

	ys = np.zeros(Lmax+1)
	LMMatrix = np.zeros((k,k))
	for j,a in enumerate(a1):
		l,m = kToLM[j]
		LMMatrix[int(l),int(m+l)] = a
		ys[int(l)]  = np.linalg.norm(LMMatrix[int(l)], ord=2)

	r_approx = Y.dot(a1)

	return a1, ys, r_approx

def simple_transform(sphericalCoordinates, Lmax, vertexAreas):
	phi, theta, radii = sphericalCoordinates
	
	for i,r,a in enumerate(zip(radii, vertexAreas)):
		for l in range(Lmax + 1):
			for m in range(-l,l+1):
				c = a*r*sph_harm(m, l, theta[i], phi[i]).real

def matlab_approach(sphericalCoordinates, Lmax, vertexAreas):
	phi, theta, radii = sphericalCoordinates
	#radii = np.ones(theta.shape[0])
	Y_N = np.ones((phi.shape[0],(Lmax+1)**2),dtype=np.complex)
	for l in range(Lmax + 1):
		Y_N[:,(l)**2:l**2+l*2+1] = np.array([sph_harm(m-l, l, theta, phi) for m in range(0, 2*l+1)]).T
	
	b = np.linalg.inv(Y_N.T.dot(Y_N)).dot(Y_N.T).dot(radii)
	return b

def processSphere(filePath):
	Lmax = int(sys.argv[2])
	vertices, faces, normals  = loadOBJ(filePath)

	bounds = Mesh(vertices.T, faces, normals).getBounds()
	p0 = [bounds[0][0],bounds[1][0],bounds[2][0],150]
	centerPoint, radius = fitSphere(vertices, p0, 150, bounds)

	sphericalCoordinates = getSphericalCoordinates(vertices, centerPoint)
	print 'Calculating areas'
	vertexAreas = np.ones(faces.shape[0])#getVoronoiArea(vertices, faces)
	vertexAreas = vertexAreas/np.linalg.norm(vertexAreas)

	phi, theta, r = sphericalCoordinates
	phi = [degrees(p) for p in phi]
	theta = [degrees(t) for t in theta]
	vertexAreas = vertexAreas / np.max(vertexAreas)

	finalYs = matlab_approach(sphericalCoordinates, Lmax, vertexAreas)

	'''
	#r = r*vertexAreas
	cirf, chi = SHExpandLSQ(r, phi, theta, Lmax)
	
	a = SHPowerSpectrum(cirf)
	#a = [np.linalg.norm(ms, ord=2) for ms in np.sum(cirf, axis=0)]
	print a
	finalYs = a

	r = r*(vertexAreas)
	sphericalCoordinates = [phi, theta, r]
	noPasses = 1

	finalYs = np.zeros(Lmax + 1)
	for i in range(noPasses):
		print "Pass " + str(i)
		a, ys, A = sirf(sphericalCoordinates, Lmax, filePath)
		finalYs += ys
		phi, theta, r = sphericalCoordinates

		r = r - A.dot(a)
		rmse = np.linalg.norm(r) / np.sqrt(len(r))
		print "RMSE: " + str(rmse)
		sphericalCoordinates = np.array([phi, theta, r])
	'''

		
	return finalYs

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print "Please specify OBJPATH NOFREQUENCIES"
		sys.exit(0)


	if sys.argv[1].endswith('.obj'):
		folderPath = ''
		files = [sys.argv[1]]
	else:
		folderPath = sys.argv[1]
		files = os.listdir(folderPath)

	numColors = 0
	for fileName in files:
		if fileName.endswith('.obj'):
			numColors += 1

	cm = plt.get_cmap('gist_rainbow')

	fig = plt.figure()
	ax = fig.add_subplot(111)
	#colors = [cm(1.*i/numColors) for i in range(numColors)]
	cNorm  = mcolors.Normalize(vmin=0, vmax=numColors-1)
	scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)

	colors = [scalarMap.to_rgba(i) for i in range(numColors)]
	lineStyles = ['-', '--', ':', '-.']
	lineWidths = [2,3]

	paddedLineStyles = np.tile(lineStyles, len(colors)/len(lineStyles))
	paddedLineWidths = np.tile(lineWidths, len(colors)/len(lineWidths))
	if len(paddedLineStyles) < len(colors):
		paddedLineStyles = np.hstack((paddedLineStyles,lineStyles[0:len(colors)-len(paddedLineStyles)]))
	if len(paddedLineWidths) < len(colors):
		paddedLineWidths = np.hstack((paddedLineWidths,lineWidths[0:len(colors)-len(paddedLineWidths)]))

	ax.set_prop_cycle(cycler('color', colors) + cycler('linestyle', paddedLineStyles) + cycler('linewidth', paddedLineWidths))
	for fileName in files:
		if fileName.endswith('.obj'):

			start_time = timeit.default_timer()
			ys = processSphere(os.path.join(folderPath,fileName))
			print(timeit.default_timer() - start_time)

			xa = np.arange(0, len(ys))
			plt.plot(xa, np.absolute(ys), label=fileName[0:-4])

			#rmse = np.linalg.norm(r - r_approx) / np.sqrt(len(a1))

	plt.legend()
	plt.grid(True)
	#plt.gca().set_yscale('log')
	plt.show()
