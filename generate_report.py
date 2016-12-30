#! /usr/bin/python

import os
import sys
import time
import importlib
import opendr
import numpy as np
from src.reporting import writeReport
from src.measure import Measure, getMeasures
from src.utils import scanLineFit
from src.shapes import Sphere
from src.fitting import distance
from src.calibration import getStereoCamerasFromCalibration, StereoCamera, cameraDistance

if __name__ == '__main__':
	
	if len(sys.argv) < 3:
		print 'Please provide as an argument the directory where the OBJ-files are located and the radius to be used'
		sys.exit(1)

	shapes = []

	for fileName in os.listdir(sys.argv[1]):
		if fileName.endswith('.obj'):
			filePath = sys.argv[1] + '/' + fileName
			print 'Loading mesh ' + fileName
			shapes.append(Sphere(filePath, float(sys.argv[2])))
	
	if len(shapes) == 0:
		print 'No meshes found!'
		sys.exit(0)

	
	measures = getMeasures(Sphere)

	cameraFocusCalibrationPath = 'calibrations/20161219172120574/'
	stereoCameras = getStereoCamerasFromCalibration(cameraFocusCalibrationPath)
	cameraPosition = (stereoCameras[7].A.position + stereoCameras[7].B.position) / 2

	baseLine = distance(stereoCameras[7].A.position, stereoCameras[7].B.position)

	#measures.append(Measure('Distance from cam', lambda x: distance(x.centerPoint,cameraPosition)))
	measures.append(Measure('Distance from cam', lambda x: cameraDistance(stereoCameras[7].A.position, stereoCameras[7].B.position, x.centerPoint)))
	#measures.append(Measure('Theoretical depth error', lambda x: (distance(x.centerPoint,cameraPosition))**2 / (baseLine* stereoCameras[7].A.focalLength)))

	results = []
	for i, shape in enumerate(shapes):
		result = []
		for j,measurement in enumerate(measures):
			result.append(measurement.execute(shape))
		results.append(result)

	
	writeReport('reportTest.xlsx', measures, results)
