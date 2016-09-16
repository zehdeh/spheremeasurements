#! /usr/bin/python

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from src.shapes import Sphere
if __name__ == '__main__':
	
	if len(sys.argv) < 2:
		print 'Please provide as an argument the directory where the OBJ-files are located'
		sys.exit(1)

	shapes = []
	for fileName in os.listdir(sys.argv[1]):
		if fileName.endswith('.obj'):
			filePath = sys.argv[1] + '/' + fileName
			print 'Loading mesh ' + fileName
			shapes.append(Sphere(filePath, 150))

	#for s in shapes:
	#	print s.totalFittingError() / len(s.vertices.T)
	

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	for s in shapes:
		ax.scatter(s.centerPoint[0], s.centerPoint[1], s.centerPoint[2])
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')
	plt.show()
