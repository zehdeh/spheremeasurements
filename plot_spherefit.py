#! /usr/bin/python

import sys
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from src.shapes import Sphere

if __name__ == '__main__':
	sphere1 = Sphere('res/tangential1/LargeSphere_tangential_01.000001.obj', 150)

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')


	ax.scatter(sphere1.vertices[0],sphere1.vertices[1],sphere1.vertices[2], color='b', marker='.')

	ax.scatter(sphere1.centerPoint[0],sphere1.centerPoint[1],sphere1.centerPoint[2], color='r',s=10)
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')
	plt.show()


