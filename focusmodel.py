#! /usr/bin/python

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

import config.defaults
from src.OBJIO import loadOBJ, loadOBJviaVTK
from src.shapes import Sphere
from src.fitting import fitSphere, fittingErrorSphere, distance, calculateMeanCurvature
from src.calibration import getStereoCamerasFromCalibration, StereoCamera, cameraDistance
from scipy.optimize import leastsq
from src.grid import getVoxelCenter

def buildModel(nominalRadius):
	cameraFocusFolderPath = 'res/final/camerafocus/cleanedup/'
	cameraFocusCalibrationPath = 'calibrations/20161219172120574/'

	stereoCameras = getStereoCamerasFromCalibration(cameraFocusCalibrationPath)
	cameraPosition = (stereoCameras[7].A.position + stereoCameras[7].B.position)/2

	baseLine = distance(stereoCameras[7].A.position, stereoCameras[7].B.position)

	i = 0
	xModel = []
	fittingErrors = []
	yModelCurv = []
	x0 = 0
	for fileName in os.listdir(cameraFocusFolderPath):
		if fileName.endswith('.obj'):
			filePath = os.path.join(cameraFocusFolderPath, fileName)

			sphere = Sphere(filePath, nominalRadius, False, loadvtk=True)
			fittingErrors.append(np.sum(np.abs(sphere.relativeFittingError())))
			yModelCurv.append(sphere.curvature.std())

			fileName = fileName[:-11]
			xModel.append(float(distance(cameraPosition, sphere.centerPoint)))

			distanceStr = fileName[fileName.find('_')+1:]
			if distanceStr == '0':
				x0 = len(xModel)-1
	
	sortedIndices = np.argsort(xModel)
	x0 = int(np.where(sortedIndices == x0)[0])
	xModel = np.asarray(xModel)[sortedIndices]
	yModelCurv = np.asarray(yModelCurv)[sortedIndices]

	return xModel, yModelCurv, x0

def fitModel(xModel, yModelCurv, x0Model, xMeasured, yMeasured):
	nModel = yModelCurv.shape[0]
	nMeasured = yMeasuredCurv.shape[0]

	fitFunction = lambda params, indices: np.abs(yMeasuredCurv[indices] - yModelCurv*params[0] - abs(params[1]))

	numStepSizesToCheck = 20

	lowestRMSE = np.inf
	maxStepSize = float(nMeasured-1) / (nModel-1)
	minStepSize = float((nMeasured)-1) / (nModel-1)*0.6
	for stepSize in np.linspace(minStepSize, maxStepSize, num=numStepSizesToCheck):
		indices = np.linspace(0, (nModel-1)*stepSize, nModel).astype(int)
		for xPos in range(0, nMeasured - indices[-1]):
			res = np.abs(yMeasuredCurv[indices + xPos] - yModelCurv)
			rmse = np.sqrt(np.mean(res**2))
			if rmse < lowestRMSE:
				lowestRMSE = rmse
				xPosLowest = xPos
				stepSizeLowest = stepSize

	indices = np.linspace(0, (nModel-1)*stepSizeLowest, nModel).astype(int) + xPosLowest
	params,flag = leastsq(fitFunction, [0,0], args=(indices,))
	yScale, yShift = params

	finalRMSE = np.sqrt(np.mean(fitFunction([yScale, yShift],indices)**2))

	return xPosLowest, stepSizeLowest, yScale, yShift

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Generates a model based on our camera focus experiment')
	parser.add_argument("--radius", help="The nominal radius to be used", type=float)
	parser.add_argument("--show-plots", help="Show plots of curves during model generation", action='store_true')
	parser.add_argument("--verbose", help="Show debug information", action='store_true')
	args = parser.parse_args()

	if args.radius is None:
		args.radius = config.defaults.nominalRadius

	if args.verbose:
		print 'Using radius ' + str(args.radius)
		print 'Building model..'

	xModel, yModelCurv, x0Model = buildModel(args.radius)
	nModel = yModelCurv.shape[0]

	if args.show_plots:
		fig = plt.figure(facecolor='white')
		plt.plot(xModel, yModelCurv)
		plt.plot([xModel[x0Model], xModel[x0Model]],[0,np.max(yModelCurv)], linestyle='dashed')
		plt.show()
		plt.xlabel('x')
		plt.title('Model')

	cameraFolderName = 'res/final/onecam_fullvolume/cleanedup'

	cameraFocusCalibrationPath = 'calibrations/20161222142605764/'
	stereoCameras = getStereoCamerasFromCalibration(cameraFocusCalibrationPath)

	gridSize = config.defaults.gridSize
	gridScale = config.defaults.gridScale
	totalCameraFocusModel = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.float)
	totalVisibilityMatrix = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.int)

	for cameraNo in stereoCameras:
		stereoCamera = stereoCameras[cameraNo]

		cameraPosition = (stereoCamera.A.position + stereoCamera.B.position) / 2

		cameraFolderName = 'res/final/onecam_fullvolume/' + str(cameraNo).zfill(2) + '/cleanedup/'
		if not os.path.exists(cameraFolderName):
			if args.verbose:
				print cameraFolderName + ' doesn\'t exist! Skipping camera ' + str(cameraNo)
			continue

		if args.verbose:
			print 'Continuing with camera ' + str(cameraNo)

		fittingErrors = []
		yMeasuredCurv = []
		xMeasured = []
		files = os.listdir(cameraFolderName)

		for fileName in files:
			if fileName.endswith('.obj'):
				filePath = os.path.join(cameraFolderName, fileName)
				sphere = Sphere(filePath, args.radius, False, loadvtk=True)

				fittingErrors.append(np.sum(np.abs(sphere.relativeFittingError())))
				yMeasuredCurv.append(sphere.curvature.std())
				
				xMeasured.append(float(distance(sphere.centerPoint, cameraPosition)))

		sortedIndices = np.argsort(xMeasured)
		xMeasured = np.array(xMeasured)[sortedIndices]
		fittingErrors = np.array(fittingErrors)[sortedIndices]
		yMeasuredCurv = np.array(yMeasuredCurv)[sortedIndices]
		nMeasured = yMeasuredCurv.shape[0]

		if args.verbose:
			iExtremeCurvatures = np.argsort(yMeasuredCurv)[-11:-1]
			print 'Models with extreme values:'
			for i in iExtremeCurvatures:
				print files[i] + ' ' + str(yMeasuredCurv[i])

		N = 10
		yMeasuredCurv = np.convolve(yMeasuredCurv, np.ones(N)/N, mode='same')
		if args.verbose:
			print 'Fitting model..'
		xPos, xStepSize, yScale, yShift = fitModel(xModel, yModelCurv, x0Model, xMeasured, yMeasuredCurv)

		yModelAdapted = yModelCurv*yScale+yShift

		indices = np.linspace(0, (nModel-1)*xStepSize, nModel).astype(int) + xPos
		xModelAbsolute = xMeasured[indices]

		if args.show_plots:
			fig = plt.figure(facecolor='white')
			plt.plot(xMeasured, yMeasuredCurv)
			plt.scatter(xModelAbsolute, yMeasuredCurv[indices], s=20,marker='o', facecolors='none', edgecolors='b')
			plt.plot(xModelAbsolute, yModelAdapted, linewidth='2')
			plt.plot([xModelAbsolute[x0Model], xModelAbsolute[x0Model]], [0,np.max(yMeasuredCurv)], linestyle='dashed')
			plt.title('Camera ' + str(cameraNo))
			plt.xlabel('x')
			plt.ylabel('std. of curvature')
			plt.ylim(0, 0.1)
			plt.show()

		visibilityMatrixFileName = os.path.join(cameraFocusCalibrationPath,\
		'visibility_' + str(cameraNo).zfill(2))
		if os.path.isfile(visibilityMatrixFileName + '.npy'):
			stereoCamera.visibilityMatrix = np.load(visibilityMatrixFileName + '.npy')
		else:
			if args.verbose:
				print 'Generating visibility matrix..'
			stereoCamera.generateVisibilityMatrix(gridSize, gridScale)
			np.save(visibilityMatrixFileName, stereoCamera.visibilityMatrix)
		totalVisibilityMatrix += stereoCamera.visibilityMatrix

		if args.verbose:
			print 'Generating model matrix..'
		cameraFocusModel = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.float)
		for i in range(gridSize[0]):
			for j in range(gridSize[1]):
				for k in range(gridSize[2]):
					if stereoCamera.visibilityMatrix[i,j,k] == 1:
						x,y,z = getVoxelCenter(gridSize, gridScale, i,j,k)

						cellDistance = distance([x,y,z], cameraPosition)
						iClosest = np.abs(xModelAbsolute-cellDistance).argmin()

						lowerBound = max(0,iClosest if cellDistance > xModelAbsolute[iClosest] else iClosest - 1)
						upperBound = min(nModel-1,iClosest if cellDistance < xModelAbsolute[iClosest] else iClosest + 1)

						cameraFocusModel[i,j,k] = (yModelAdapted[lowerBound]+yModelAdapted[upperBound])/2
		totalCameraFocusModel += cameraFocusModel
	np.save(os.path.join(cameraFocusCalibrationPath,'total_visibility'), totalVisibilityMatrix)

	totalCameraFocusModel[np.nonzero(totalVisibilityMatrix)] /= totalVisibilityMatrix[np.nonzero(totalVisibilityMatrix)]

	totalVisibilityMatrix = totalVisibilityMatrix.astype(np.float)/totalVisibilityMatrix.max()
	totalCameraFocusModel *= totalVisibilityMatrix

	np.save(os.path.join(cameraFocusCalibrationPath,'focusmodel'),totalCameraFocusModel)
