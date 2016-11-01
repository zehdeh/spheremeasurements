#! /usr/bin/python

import os
import sys
import timeit
import numpy as np
import time

from math import isnan, sin, degrees
from scipy.special import sph_harm
from src.OBJIO import loadOBJ
from src.fitting import fitSphere
import src.spherical_harmonics as sh
from src.mesh import Mesh
from src.voronoi import getVoronoiArea
from scipy.optimize import leastsq
from scipy.linalg import lstsq
from pymesh import form_mesh

from cycler import cycler
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as mplcm

def simple_transform(sphericalCoordinates, Lmax, vertexAreas, removeCoefficients=False):
	phi, theta, radii = sphericalCoordinates
	
	totalArea = np.sum(vertexAreas)
	ys = np.zeros(Lmax + 1)
	coefficients = np.zeros((Lmax+1)**2, dtype=np.complex)
	for l in range(Lmax + 1):
		for m in range(-l,l+1):
			j = l**2 + l + m
			coefficients[j] = np.sum(radii*sph_harm(m, l, phi, theta)*vertexAreas)/totalArea
			if removeCoefficients:
				radii -= np.absolute(coefficients[j]*sph_harm(m, l, phi, theta))
		ys[l] = np.linalg.norm(coefficients[l**2:(l**2 + 2*l + 1)], ord=2)
	
	return ys

def matlab_approach(sphericalCoordinates, Lmax, vertexAreas):
	phi, theta, radii = sphericalCoordinates
	#radii = np.ones(theta.shape[0])
	Y_N = np.ones((phi.shape[0],(Lmax+1)**2),dtype=np.complex)
	for l in range(Lmax + 1):
		Y_N[:,(l)**2:l**2+l*2+1] = np.array([sph_harm(m-l, l, phi, theta) for m in range(0, 2*l+1)]).T
	
	#b = np.linalg.inv(Y_N.T.dot(Y_N)).dot(Y_N.T).dot(radii)
	b = lstsq(Y_N*np.sqrt(vertexAreas)[:,None],radii*np.sqrt(vertexAreas))[0]
	return Y_N,b

def processSphere(filePath):
	print 'Processing ' + os.path.split(filePath)[-1]
	Lmax = int(sys.argv[2])
	vertices, faces, normals  = loadOBJ(filePath)

	bounds = Mesh(vertices.T, faces, normals).getBounds()
	p0 = [bounds[0][0],bounds[1][0],bounds[2][0],150]
	centerPoint, radius = fitSphere(vertices, p0, 150, bounds)

	sphericalCoordinates = sh.getSphericalCoordinates(vertices, centerPoint)
	print 'Calculating areas'
	mesh = form_mesh(vertices, faces)
	mesh.add_attribute('vertex_area')
	vertexAreas = mesh.get_attribute('vertex_area')
	vertexAreas = vertexAreas/np.linalg.norm(vertexAreas)

	meanArea = np.mean(vertexAreas)
	stdArea = np.std(vertexAreas)
	
	excludedAreasIdx = vertexAreas > (meanArea + stdArea)
	sphericalCoordinates = sphericalCoordinates[:,~excludedAreasIdx]
	vertexAreas = vertexAreas[~excludedAreasIdx]

	phi, theta, r = sphericalCoordinates
	phi = [degrees(p) for p in phi]
	theta = [degrees(t) for t in theta]

	vertexAreas = vertexAreas / np.linalg.norm(vertexAreas, ord=2)
	
	#finalYs = simple_transform(sphericalCoordinates, Lmax, vertexAreas, removeCoefficients=True)

	noPasses = 1
	finalA = np.zeros((Lmax+1)**2, dtype=np.complex)
	for i in range(noPasses):
		A,a = matlab_approach(sphericalCoordinates, Lmax, vertexAreas)
		finalA += a

		r = r - A.dot(a)
		rmse = np.linalg.norm(r) / np.sqrt(len(r))
		print "RMSE: " + str(rmse)
		r = np.absolute(r)
		sphericalCoordinates = np.array([phi, theta, r])

	#r = r*vertexAreas

	'''
	r = r*(vertexAreas)
	sphericalCoordinates = [phi, theta, r]

	finalYs = np.zeros(Lmax + 1)
		print "Pass " + str(i)
		a, ys, A = sh.sirf(sphericalCoordinates, Lmax, filePath, vertexAreas)
		finalYs += ys
		finalYs = a
		phi, theta, r = sphericalCoordinates

		r = r - A.dot(a)
		rmse = np.linalg.norm(r) / np.sqrt(len(r))
		print "RMSE: " + str(rmse)
		sphericalCoordinates = np.array([phi, theta, r])

		

	'''
	finalYs = np.zeros(Lmax+1)
	for i in range(Lmax+1):
		lowerBound = i**2
		upperBound = (i+1)**2
		finalYs[i] = np.linalg.norm(finalA[lowerBound:upperBound], ord=2)
	
	print 'Summed coefficients: ' + str(np.sum(finalYs))

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

	plt.legend(fontsize=9)
	plt.grid(True)
	#plt.gca().set_yscale('log')
	plt.show()
