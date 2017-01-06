#! /usr/bin/python

import sys
import os
import numpy as np
import argparse

import config.defaults
from src.OBJIO import loadOBJ
from src.fitting import fitSphere, fittingErrorSphere, distance

def checkDir(directory):
	if not os.path.isdir(directory):
		raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(directory))
	return directory

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Generates a spreadsheet')
	parser.add_argument("folder", help="The folder with the OBJ files", type=checkDir)
	parser.add_argument("--verbose", help="Show debug information", action='store_true')
	parser.add_argument("--dry-run", help="Do not actually delete anything", action='store_true')

	args = parser.parse_args()

	nominalRadius = config.defaults.nominalRadius
	
	nEliminated = 0
	nTotal = 0

	for fileName in os.listdir(args.folder):
		if fileName.endswith('.obj'):
			nTotal += 1

			filePath = os.path.join(args.folder, fileName)
			vertices, faces, normals = loadOBJ(filePath)

			deleting = False
			if vertices.shape[0] < config.defaults.minNumberVertices:
				if args.verbose:
					print 'Vertices: ' + str(vertices.shape[0])
					print 'Deleting ' + filePath
				nEliminated += 1
				if not args.dry_run:
					os.remove(filePath)
				continue

			centerPoint, fittedRadius = fitSphere(vertices, nominalRadius, False)
			fittingError = np.sum(np.fabs(nominalRadius - distance(vertices.T, centerPoint)))/vertices.shape[0]

			centerPoint, fittedRadius = fitSphere(vertices, nominalRadius, True, True)

			if fittingError > 5.5:
				if args.verbose:
					print 'Fitting error: ' + str(fittingError)
					print 'Deleting ' + filePath

				deleting = True
			if (fittedRadius/nominalRadius) > 1.2 or (fittedRadius/nominalRadius) < 0.8:
				if args.verbose:
					print 'Fitted Radius: ' + str(fittedRadius)
					print 'Deleting ' + filePath
				deleting = True

			if deleting:
				nEliminated += 1
				if not args.dry_run:
					os.remove(filePath)
	
	print 'Deleted ' + str(nEliminated) + ' files in total.'
	print str(nTotal - nEliminated) + ' files left.'
