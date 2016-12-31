#! /usr/bin/python

import sys
import os
import numpy as np
from src.OBJIO import loadOBJ
from src.fitting import fitSphere, fittingErrorSphere, distance

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print 'Please provide a folder name and the nominal radius.'
		sys.exit(0)

	folderName = sys.argv[1]
	nominalRadius = float(sys.argv[2])

	for fileName in os.listdir(folderName):
		if fileName.endswith('.obj'):
			filePath = os.path.join(folderName, fileName)
			vertices, faces, normals = loadOBJ(filePath)

			centerPoint, fittedRadius = fitSphere(vertices, nominalRadius, False)
			fittingError = np.sum(np.fabs(nominalRadius - distance(vertices.T, centerPoint)))/vertices.shape[0]

			centerPoint, fittedRadius = fitSphere(vertices, nominalRadius, True)

			if vertices.shape[0] < 250 or fittingError > 5.5 or (fittedRadius/nominalRadius) > 1.1 or (fittedRadius/nominalRadius) < 0.9:
				print 'Fitted Radius: ' + str(fittedRadius)
				print 'Vertices: ' + str(vertices.shape[0])
				print 'Fitting error: ' + str(fittingError)
				print 'Deleting ' + filePath

				os.remove(filePath)
