#! /usr/bin/python

import sys
import os
import argparse
import numpy as np

from PyQt5.QtWidgets import QApplication
from src.ui import MainWindow
from src.calibration import getStereoCamerasFromCalibration, StereoCamera
from scipy.ndimage.filters import convolve
from src.shapes import Sphere
from src.grid import generateCurvatureGrid, generateVectorField, generateErrorGrid

def checkDir(directory):
	if not os.path.isdir(directory):
		raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(directory))
	return directory

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Displays the error as a volume rendering')
	parser.add_argument("folder", help="The folder with the OBJ files", type=checkDir)
	parser.add_argument("calibration", help="If you are interested in camera distance, provide a calibration folder", type=checkDir)
	parser.add_argument("--verbose", help="Show debug information", action='store_true')
	parser.add_argument("--rebuild-cache", help="Rebuild the cache files", action='store_true')
	args = parser.parse_args()

	app = QApplication(sys.argv)

	gridSize = [50,50,50]
	scannerVolumeSize = [2400,3000,2400]
	gridScale = [scannerVolumeSize[0] / gridSize[0], scannerVolumeSize[1] / gridSize[1], scannerVolumeSize[2] / gridSize[2]]

	vectorCacheFileName = os.path.join(args.folder,'vectorField' + '_' + str(gridSize[0]) + '_' + str(gridSize[1]) + '_' + str(gridSize[2]))
	errorCacheFileName = os.path.join(args.folder,'fittingError' + '_' + str(gridSize[0]) + '_' + str(gridSize[1]) + '_' + str(gridSize[2]))
	nMatrixCacheFileName = os.path.join(args.folder,'nMatrix' + '_' + str(gridSize[0]) + '_' + str(gridSize[1]) + '_' + str(gridSize[2]))

	if os.path.isfile(errorCacheFileName + '.npy') and not args.rebuild_cache:
		if args.verbose:
			print 'Using cache..'
		nMatrix = np.load(nMatrixCacheFileName + '.npy')
		errorMatrix = np.load(errorCacheFileName + '.npy')
		vectorField = np.load(vectorCacheFileName + '.npy')
	else:
		spheres = []
		for fileName in os.listdir(args.folder):
			if fileName.endswith('.obj'):
				if args.verbose:
					print 'Loading ' + fileName[0:-4]
				nominalRadius = 80.06505
				spheres.append(Sphere(os.path.join(args.folder, fileName), nominalRadius))

		if len(spheres) == 0:
			raise RuntimeError('No scans were found')

		if args.verbose:
			print 'Generating error matrix...'
		errorMatrix, nMatrix = generateErrorGrid(gridSize, gridScale, spheres, args.verbose)
		if args.verbose:
			print 'Generating vector field...'
		vectorField = generateVectorField(gridSize, gridScale, spheres, nMatrix, errorMatrix)

		np.save(errorCacheFileName, errorMatrix)
		np.save(vectorCacheFileName, vectorField)
		np.save(nMatrixCacheFileName, nMatrix)

	confidenceMatrix = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.float)
	confidenceMatrix[np.nonzero(nMatrix)] = nMatrix[np.nonzero(nMatrix)]

	confidenceMatrix = confidenceMatrix/confidenceMatrix.max()

	errorMatrix[np.nonzero(nMatrix)] = errorMatrix[np.nonzero(nMatrix)] / errorMatrix.max()
	vectorField *= 10

	calibrationFolderPath = sys.argv[2]
	stereoCameras = getStereoCamerasFromCalibration(calibrationFolderPath)
	
	window = MainWindow(stereoCameras, gridSize, gridScale, errorMatrix, vectorField, confidenceMatrix, args.verbose)
	sys.exit(app.exec_())
