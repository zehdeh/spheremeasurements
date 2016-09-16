#! /usr/bin/python

import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from src.shapes import Sphere

if __name__ == '__main__':
	sphere = Sphere('res/tangential1/LargeSphere_tangential_01.000001.obj', 150)
	sphere2 = Sphere('res/tangential1/LargeSphere_tangential_01.000016.obj', 150)
	sphere3 = Sphere('res/tangential1/LargeSphere_tangential_01.000040.obj', 150)

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	ax.scatter(sphere.vertices[0],sphere.vertices[1],sphere.vertices[2], color='b', marker='.')
	ax.scatter(sphere2.vertices[0],sphere2.vertices[1],sphere2.vertices[2], color='g', marker='.')
	ax.scatter(sphere3.vertices[0],sphere3.vertices[1],sphere3.vertices[2], color='b', marker='.')

	ax.scatter(sphere.centerPoint[0],sphere.centerPoint[1],sphere.centerPoint[2], color='r',s=10)
	ax.scatter(sphere2.centerPoint[0],sphere2.centerPoint[1],sphere2.centerPoint[2], color='r',s=10)
	ax.scatter(sphere3.centerPoint[0],sphere3.centerPoint[1],sphere3.centerPoint[2], color='r',s=10)
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')
	ax.scatter(sphere.centerPoint[0], sphere.centerPoint[1], sphere.centerPoint[2], color='r')
	plt.show()


