#! /usr/bin/python

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from src.shapes import Sphere
from src.utils import scanLineFit, scanCircleFit


if __name__ == '__main__':
	
	if len(sys.argv) < 2:
		print 'Please provide as an argument the directory where the OBJ-files are located'
		sys.exit(1)

	files = os.listdir(sys.argv[1])

	shapes = []
	for fileName in files:
		if fileName.endswith('.obj'):
			filePath = sys.argv[1] + '/' + fileName
			print 'Loading mesh ' + fileName
			shapes.append(Sphere(filePath, 150))

	centerPoints = np.asarray([s.centerPoint for s in shapes])
	#for s in shapes:
	#	print s.totalFittingError() / len(s.vertices.T)

	fittingErrors = np.asarray([s.totalFittingError() / s.vertices.shape[1] for s in shapes])
	indices = fittingErrors.argsort(0)
	print type(indices)

	for f,e in zip(np.asarray(files)[indices], fittingErrors[indices]):
		print str(f) + ' ' + str(e)


	#mean, unit_v, projectedPoints, pointDistances = scanLineFit(centerPoints)
	p1 = scanCircleFit(centerPoints)

	#linePoints = unit_v * np.mgrid[-800:800:2j][:,np.newaxis]
	#linePoints += mean

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	centerPoints = centerPoints.T
	ax.scatter(centerPoints[0], centerPoints[1], centerPoints[2])
	#ax.scatter(projectedPoints[0], projectedPoints[1], projectedPoints[2], c='r', marker='x', s=100)
	ax.scatter(p1[0], p1[1], p1[2], c='r', marker='x', s=100)
	#ax.plot3D(*linePoints.T)
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')

	# Create cubic bounding box to simulate equal aspect ratio
	max_range = np.array([centerPoints[0].max()-centerPoints[0].min(), centerPoints[1].max()-centerPoints[1].min(), centerPoints[2].max()-centerPoints[2].min()]).max()
	Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(centerPoints[0].max()+centerPoints[0].min())
	Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(centerPoints[1].max()+centerPoints[1].min())
	Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(centerPoints[2].max()+centerPoints[2].min())
	# Comment or uncomment following both lines to test the fake bounding box:
	for xb, yb, zb in zip(Xb, Yb, Zb):
		ax.plot([xb], [yb], [zb], 'w')
	
	plt.show()
