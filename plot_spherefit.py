#! /usr/bin/python

import sys
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from src.shapes import Sphere

if __name__ == '__main__':
	sphere1 = Sphere('res/tangential1/LargeSphere_tangential_01.000001.obj', 150)
	sphere2 = Sphere('res/tangential1/LargeSphere_tangential_01.000025.obj', 150)
	sphere3 = Sphere('res/tangential1/LargeSphere_tangential_01.000049.obj', 150)

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')


	ax.scatter(sphere1.vertices[0],sphere1.vertices[1],sphere1.vertices[2], color='b', marker='.')
	ax.scatter(sphere2.vertices[0],sphere2.vertices[1],sphere2.vertices[2], color='b', marker='.')
	ax.scatter(sphere3.vertices[0],sphere3.vertices[1],sphere3.vertices[2], color='b', marker='.')

	ax.scatter(sphere1.centerPoint[0],sphere1.centerPoint[1],sphere1.centerPoint[2], color='r',s=10)
	ax.scatter(sphere2.centerPoint[0],sphere2.centerPoint[1],sphere2.centerPoint[2], color='r',s=10)
	ax.scatter(sphere3.centerPoint[0],sphere3.centerPoint[1],sphere3.centerPoint[2], color='r',s=10)
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')
	plt.show()


