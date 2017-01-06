#! /usr/bin/python

import sys
import os
import numpy as np

import config.defaults
from src.OBJIO import loadOBJ
from src.fitting import fitSphere, fittingErrorSphere, distance

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print 'Please provide a folder'
		sys.exit(0)

	folderName = sys.argv[1]
	nominalRadius = config.defaults.nominalRadius
	
	nEliminated = 0
	nTotal = 0

	for fileName in os.listdir(folderName):
		if fileName.endswith('.obj'):
			nTotal += 1
			filePath = os.path.join(folderName, fileName)
			vertices, faces, normals = loadOBJ(filePath)
			if vertices.shape[0] < config.defaults.minNumberVertices:
				print 'Vertices: ' + str(vertices.shape[0])
				print 'Deleting ' + filePath
				os.remove(filePath)
				nEliminated += 1
				continue

			centerPoint, fittedRadius = fitSphere(vertices, nominalRadius, False)
			fittingError = np.sum(np.fabs(nominalRadius - distance(vertices.T, centerPoint)))/vertices.shape[0]

			centerPoint, fittedRadius = fitSphere(vertices, nominalRadius, True, True)

			if fittingError > 5.5 or (fittedRadius/nominalRadius) > 1.6 or (fittedRadius/nominalRadius) < 0.6:
				print 'Fitted Radius: ' + str(fittedRadius)
				print 'Fitting error: ' + str(fittingError)
				print 'Deleting ' + filePath

				nEliminated += 1

				os.remove(filePath)
	
	print 'Deleted ' + str(nEliminated) + ' files in total.'
	print str(nTotal - nEliminated) + ' files left.'
