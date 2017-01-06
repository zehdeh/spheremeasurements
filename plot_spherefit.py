#! /usr/bin/python

import sys
import matplotlib.pyplot as plt
import numpy as np

import config.defaults
from mpl_toolkits.mplot3d import Axes3D
from src.OBJIO import loadOBJviaVTK
from src.fitting import fitSphere, distance
from src.shapes import Sphere

if __name__ == '__main__':

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	nominalRadius = config.defaults.nominalRadius
	sphere = Sphere(sys.argv[1], nominalRadius, False)

	ax.scatter(sphere.vertices.T[0],sphere.vertices.T[1],sphere.vertices.T[2], color='b', marker='.')

	print 'Position of centerPoint: ' + str(sphere.centerPoint)
	print 'Fitted radius: ' + str(sphere.fittedRadius)
	print 'Relative fitting error: ' + str(np.mean(sphere.relativeFittingError()))

	errors = np.abs(sphere.fittingError())
	iLargestError = np.argmax(errors)
	largestError = sphere.vertices[iLargestError]
	print 'Point with largest error: ' + str(iLargestError) + ': ' + str(errors.max()) + ' ' + str(largestError)
	print distance(sphere.centerPoint, sphere.vertices[iLargestError])

	ax.scatter(sphere.centerPoint[0],sphere.centerPoint[1],sphere.centerPoint[2], color='r',s=10)

	ax.scatter(largestError[0],largestError[1],largestError[2], color='black',s=20)
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')
	plt.show()

