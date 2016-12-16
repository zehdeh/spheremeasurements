#! /usr/bin/python

import sys
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from src.vertexarea import getVertexAreas

from src.OBJIO import loadOBJ, getBounds
from src.fitting import fitSphere
import src.spherical_harmonics as sh
import numpy as np

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print 'Please specify OBJ file and LMax'
		sys.exit(0)

	fileName = sys.argv[1]
	if not os.path.isfile(fileName):
		print 'File doesn\'t exist!'
		sys.exit(0)

	if not fileName.endswith('.obj'):
		print 'Not an OBJ file!'
		sys.exit(0)

	vertices, faces, normals = loadOBJ(fileName)
	bounds = getBounds(vertices.T)
	p0 = [bounds[0][0],bounds[1][0],bounds[2][0],150]
	centerPoint, radius = fitSphere(vertices, p0, 150, bounds)

	sphericalCoordinates = sh.getSphericalCoordinates(vertices, centerPoint)
	vertexAreas = getVertexAreas(faces, vertices)

	Lmax = int(sys.argv[2])

	phi, theta, radii = sphericalCoordinates
	finalYs, coefficients = sh.simple_transform(sphericalCoordinates, Lmax, vertexAreas)

	reconstructedRadii = sh.back_transform(sphericalCoordinates, coefficients, Lmax, vertexAreas)

	coords = sh.getCartesianCoordinates(sphericalCoordinates[0], sphericalCoordinates[1], reconstructedRadii, centerPoint)

	fig = plt.figure()

	ax = fig.add_subplot(111, projection='3d')

	ax.scatter(coords.T[0],coords.T[1],coords.T[2], color='b', marker='.')

	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')
	plt.show()

	print finalYs
	print reconstructedRadii
	print radii
	r = radii - reconstructedRadii
	rmse = np.linalg.norm(r) / np.sqrt(r.shape[0])
	print "RMSE: " + str(rmse)
