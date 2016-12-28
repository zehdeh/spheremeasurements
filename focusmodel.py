#! /usr/bin/python

import os
import sys
from src.OBJIO import loadOBJ
from src.fitting import fitSphere, fittingErrorSphere, distance
from src.calibration import getStereoCamerasFromCalibration, StereoCamera

if __name__ == '__main__':
	cameraFocusFolderPath = 'res/final/camerafocus/cleanedup/'
	cameraFocusCalibrationPath = 'calibrations/20161219172120574/'
	nominalRadius = 80.065605

	stereoCameras = getStereoCamerasFromCalibration(cameraFocusCalibrationPath)
	cameraPosition = (stereoCameras[7].A.position + stereoCameras[7].B.position)/2

	baseLine = distance(stereoCameras[7].A.position, stereoCameras[7].B.position)
	print baseLine
	print stereoCameras[7].A.focalLength
	print stereoCameras[7].B.focalLength

	for fileName in os.listdir(cameraFocusFolderPath):
		if fileName.endswith('.obj'):
			vertices, faces, normals = loadOBJ(os.path.join(cameraFocusFolderPath, fileName))

			centerPoint, fittedRadius = fitSphere(vertices, nominalRadius)
			fittingError = fittingErrorSphere(centerPoint.tolist() + [fittedRadius], vertices)
