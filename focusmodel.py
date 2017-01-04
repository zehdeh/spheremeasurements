#! /usr/bin/python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from src.OBJIO import loadOBJ, loadOBJviaVTK
from src.shapes import Sphere
from src.fitting import fitSphere, fittingErrorSphere, distance, calculateMeanCurvature
from src.calibration import getStereoCamerasFromCalibration, StereoCamera, cameraDistance
from scipy.optimize import least_squares

def fitPolynomial(params, fittingErrors, fittedPolynomial):
	polyX = np.linspace(params[0], params[1], len(fittingErrors))
	return np.fabs(fittingErrors - fittedPolynomial(polyX)*params[2])

if __name__ == '__main__':
	cameraFocusFolderPath = 'res/final/camerafocus/cleanedup/'
	cameraFocusCalibrationPath = 'calibrations/20161219172120574/'
	nominalRadius = 80.065605

	stereoCameras = getStereoCamerasFromCalibration(cameraFocusCalibrationPath)
	cameraPosition = (stereoCameras[7].A.position + stereoCameras[7].B.position)/2

	baseLine = distance(stereoCameras[7].A.position, stereoCameras[7].B.position)

	i = 0
	dist = []
	fittingErrors = []
	curvatureStd = []
	zeroPosition = 0
	for fileName in os.listdir(cameraFocusFolderPath):
		if fileName.endswith('.obj'):
			filePath = os.path.join(cameraFocusFolderPath, fileName)

			sphere = Sphere(filePath, nominalRadius, False)
			fittingErrors.append(np.sum(np.abs(sphere.relativeFittingError())))
			curvatureStd.append(sphere.curvature.std())

			fileName = fileName[:-11]
			dist.append(float(distance(cameraPosition, sphere.centerPoint)))

			distanceStr = fileName[fileName.find('_')+1:]
			if distanceStr == '0':
				zeroPosition = dist[-1]
	
	sortedIndices = np.argsort(dist)
	dist = np.asarray(dist)[sortedIndices]
	fittingErrors = np.asarray(fittingErrors)[sortedIndices]
	curvatureStd = np.asarray(curvatureStd)[sortedIndices]

	polyCoefficients = np.polyfit(dist, fittingErrors, 3)
	fittedPolynomial = np.poly1d(polyCoefficients)

	fig = plt.figure(1)

	plt.plot(dist, fittingErrors)
	plt.plot(dist, curvatureStd)
	plt.plot([zeroPosition, zeroPosition],[0,0.15])
	plt.plot(dist, fittedPolynomial(dist), color='r')
	plt.show()
	originalScale = dist

	cameraFolderName = 'res/final/onecam_fullvolume/cleanedup'

	cameraFocusCalibrationPath = 'calibrations/20161222142605764/'
	stereoCameras = getStereoCamerasFromCalibration(cameraFocusCalibrationPath)

	gridSize = [50,50,50]
	scannerVolumeSize = [3000,3000,3000]
	gridScale = [scannerVolumeSize[0] / gridSize[0], scannerVolumeSize[1] / gridSize[1], scannerVolumeSize[2] / gridSize[2]]

	for cameraNo in stereoCameras:
		stereoCamera = stereoCameras[cameraNo]

		cameraPosition = (stereoCamera.A.position + stereoCamera.B.position) / 2

		cameraFolderName = 'res/final/onecam_fullvolume/' + str(cameraNo).zfill(2) + '/cleanedup/'
		if not os.path.exists(cameraFolderName):
			print cameraFolderName + ' doesn\'t exist! Skipping camera ' + str(cameraNo)
			continue

		print 'Continuing here!'

		fittingErrors = []
		curvatureStd = []
		dist = []

		for fileName in os.listdir(cameraFolderName):
			if fileName.endswith('.obj'):
				filePath = os.path.join(cameraFolderName, fileName)
				#vertices, faces, normals = loadOBJ(filePath)
				sphere = Sphere(filePath, nominalRadius, False)

				fittingErrors.append(np.sum(np.abs(sphere.relativeFittingError())))
				curvatureStd.append(sphere.curvature.std())
				
				dist.append(float(distance(sphere.centerPoint, cameraPosition)))

		sortedIndices = np.argsort(dist)
		dist = np.array(dist)[sortedIndices]
		fittingErrors = np.array(fittingErrors)[sortedIndices]
		curvatureStd = np.array(curvatureStd)[sortedIndices]

		params0 = [originalScale[0], originalScale[-1], 1, 0]
		res = least_squares(fitPolynomial, params0, bounds=([originalScale[0], originalScale[0], 0, 0], [originalScale[-1], originalScale[-1], 100, 100]),args=(fittingErrors, fittedPolynomial))

		lowerBound, upperBound, yScale, yShift = res.x
		polyX = np.linspace(lowerBound, upperBound, fittingErrors.shape[0])

		print dist[0]
		print dist[-1]
		print 'lb: ' + str(lowerBound)
		print 'ub: ' + str(upperBound)
		print 'estimated x0: ' + str(zeroPosition)

		scale = (dist[-1]-dist[0])/(upperBound - lowerBound)
		xAbsOffset = (zeroPosition - lowerBound)*scale
		xAbs = dist[0] + xAbsOffset
		print 'scale: ' + str(scale)
		print 'x0a: ' + str(xAbs)

		fig2 = plt.figure(2)
		#plt.plot(dist, fittingErrors)
		plt.plot(dist, curvatureStd)
		#plt.plot(dist, fittedPolynomial(polyX)*yScale + yShift)
		#plt.plot([xAbs,xAbs], [0,5])
		plt.show()

		visibilityMatrixFileName = os.path.join(cameraFocusCalibrationPath,\
		'visibility_' + str(cameraNo).zfill(2))
		if os.path.isfile(visibilityMatrixFileName + '.npy'):
			stereoCamera.visibilityMatrix = np.load(visibilityMatrixFileName + '.npy')
		else:
			print 'Generating visibility'
			stereoCamera.generateVisibilityMatrix(gridSize, gridScale)
			np.save(visibilityMatrixFileName, stereoCamera.visibilityMatrix)

		cameraFocusModel = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.float)
		for i in range(gridSize[0]):
			for j in range(gridSize[1]):
				for k in range(gridSize[2]):
					if stereoCamera.visibilityMatrix[i,j,k] == 1:
						x = (i*gridScale[0]) - (gridSize[0]*gridScale[0]/2) + gridScale[0]/2
						y = (j*gridScale[1]) - (gridSize[1]*gridScale[1]/2) + gridScale[1]/2
						z = (k*gridScale[2]) - (gridSize[2]*gridScale[2]/2) + gridScale[2]/2

						#print distance([x,y,z], cameraPosition)
		#print len(np.nonzero(stereoCamera.visibilityMatrix))
		#	print i, j, k

