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
from src.vertexarea import getVertexAreas
from scipy.optimize import leastsq
from scipy.linalg import lstsq

from cycler import cycler
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as mplcm


def processSphere(filePath):
	print 'Processing ' + os.path.split(filePath)[-1][:-4]
	Lmax = int(sys.argv[2])

	radiusNominal = 80

	cacheFileName = filePath[:-4] + '_cached_' + str(Lmax)
	if os.path.isfile(cacheFileName + '.npy'):
		finalYs = np.load(cacheFileName + '.npy')
	else:
		vertices, faces, normals  = loadOBJ(filePath)

		bounds = Mesh(vertices.T, faces, normals).getBounds()
		p0 = [bounds[0][0],bounds[1][0],bounds[2][0],radiusNominal]
		centerPoint, radius = fitSphere(vertices, p0, radiusNominal, bounds)

		sphericalCoordinates = sh.getSphericalCoordinates(vertices, centerPoint)
		print 'Calculating areas'

		vertexAreas = getVertexAreas(faces, vertices)

		print 'Finished calculating areas'

		finalYs, coefficients = sh.simple_transform(sphericalCoordinates, Lmax, vertexAreas)
		np.save(cacheFileName, finalYs)

	'''
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

	finalYs = np.zeros(Lmax+1)
	for i in range(Lmax+1):
		lowerBound = i**2
		upperBound = (i+1)**2
		finalYs[i] = np.linalg.norm(finalA[lowerBound:upperBound], ord=2)
	
	print 'Summed coefficients: ' + str(np.sum(finalYs))
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
			print 'total: ' + str(np.sum(ys))
			print(timeit.default_timer() - start_time)

			xa = np.arange(0, len(ys))
			plt.plot(xa, ys, label=fileName[0:-4])

	plt.legend(fontsize=9)
	plt.grid(True)
	#plt.gca().set_yscale('log')
	plt.show()
