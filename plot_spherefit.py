#! /usr/bin/python

import sys
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from src.OBJIO import loadOBJviaVTK
from src.fitting import fitSphere, getBounds

if __name__ == '__main__':
	vertices, faces, normals, polyData = loadOBJviaVTK(sys.argv[1])

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	radiusNominal = int(sys.argv[2])

	centerPoint, radius = fitSphere(vertices, radiusNominal)

	ax.scatter(vertices.T[0],vertices.T[1],vertices.T[2], color='b', marker='.')

	ax.scatter(centerPoint[0],centerPoint[1],centerPoint[2], color='r',s=10)
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')
	plt.show()

	print centerPoint
