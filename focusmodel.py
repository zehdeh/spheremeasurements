#! /usr/bin/python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from src.OBJIO import loadOBJ, loadOBJviaVTK
from src.fitting import fitSphere, fittingErrorSphere, distance, calculateMeanCurvature
from src.calibration import getStereoCamerasFromCalibration, StereoCamera, cameraDistance
from scipy.optimize import least_squares

def fitPoly(params, fittingErrors, fittedPoly, originalScale):
	return (fittingErrors[params[4]:len(originalScale)+params[4]] - params[0]*fittedPoly(params[1]*(originalScale-params[2])) + params[3])

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
	for fileName in os.listdir(cameraFocusFolderPath):
		if fileName.endswith('.obj'):
			filePath = os.path.join(cameraFocusFolderPath, fileName)

			vertices, faces, normals = loadOBJ(filePath)

			centerPoint, fittedRadius = fitSphere(vertices, nominalRadius, False)
			fittingErrors.append(np.sum(np.fabs(nominalRadius - distance(vertices.T, centerPoint)))/vertices.shape[0])

			fileName = fileName[:-11]
			distanceStr = fileName[fileName.find('_')+1:]
			if 'P' in distanceStr:
				dist.append(int(distanceStr[1:]))
			elif 'M' in distanceStr:
				dist.append(-int(distanceStr[1:]))
			else:
				dist.append(0)
			i += 1
	
	sortedIndices = np.argsort(dist)
	dist = np.asarray(dist)[sortedIndices]*10
	fittingErrors = np.asarray(fittingErrors)[sortedIndices]

	polyCoefficients = np.polyfit(dist, fittingErrors, 3)
	fittedPolynomial = np.poly1d(polyCoefficients)

	fig = plt.figure(1)

	plt.plot(dist, fittingErrors)
	plt.plot(dist, fittedPolynomial(dist), color='r')
	plt.show()
	originalScale = dist

	cameraFolderName = 'res/final/onecam_fullvolume/cleanedup'

	cameraFocusCalibrationPath = 'calibrations/20161222142605764/'
	stereoCameras = getStereoCamerasFromCalibration(cameraFocusCalibrationPath)

	for cameraNo in stereoCameras:
		stereoCamera = stereoCameras[cameraNo]

		cameraPosition = (stereoCamera.A.position + stereoCamera.B.position) / 2

		cameraFolderName = 'res/final/onecam_fullvolume/' + str(cameraNo).zfill(2) + '/cleanedup/'
		if not os.path.exists(cameraFolderName):
			print cameraFolderName + ' doesn\'t exist! Skipping camera ' + str(cameraNo)
			continue

		print 'Continuing here!'

		fittingErrors = []
		meanCurvatureStd = []
		dist = []

		for fileName in os.listdir(cameraFolderName):
			if fileName.endswith('.obj'):
				filePath = os.path.join(cameraFolderName, fileName)
				#vertices, faces, normals = loadOBJ(filePath)

				vertices, faces, normals, polyMesh = loadOBJviaVTK(filePath)
				meanCurvature = calculateMeanCurvature(polyMesh)
				meanCurvatureStd.append(np.std(meanCurvature))

				centerPoint, fittedRadius = fitSphere(vertices, nominalRadius, False)
				fittingErrors.append(np.sum(np.fabs(nominalRadius - distance(vertices.T, centerPoint)))/vertices.shape[0])
				
				dist.append(float(distance(centerPoint, cameraPosition)))

		sortedIndices = np.argsort(dist)
		dist = np.array(dist)[sortedIndices]
		fittingErrors = np.array(fittingErrors)[sortedIndices]
		meanCurvatureStd = np.array(meanCurvatureStd)[sortedIndices]

		polyCoefficients = np.polyfit(dist, fittingErrors, 3)
		fittedPolynomial = np.poly1d(polyCoefficients)

		fig2 = plt.figure(2)
		plt.plot(dist, fittingErrors)
		plt.plot(dist, fittedPolynomial(dist))
		plt.plot(dist, 50*meanCurvatureStd)


		plt.show()
