#! /usr/bin/python

import os
import sys
import time
import importlib
import opendr
import numpy as np
from src.reporting import writeReport
from src.measure import Measure, getMeasures
from src.utils import scanCurveFit

if __name__ == '__main__':
	
	if len(sys.argv) < 3:
		print 'Please provide as an argument the directory where the OBJ-files are located and the shape to be used'
		sys.exit(1)

	module = importlib.import_module('src.shapes')
	class_ = getattr(module, sys.argv[2])
	
	shapes = []

	for fileName in os.listdir(sys.argv[1]):
		if fileName.endswith('.obj'):
			filePath = sys.argv[1] + '/' + fileName
			print 'Loading mesh ' + fileName
			shapes.append(class_(filePath, *sys.argv[3:]))

	
	#measures = getMeasures(Sphere)
	measures = getMeasures(class_)

	centerPoints = [s.centerPoint for s in shapes]
	mean, unit_v, projectedPoints, pointDistances = scanCurveFit(centerPoints)

	measures.append(Measure('Distance from cam', lambda x: (x.centerPoint - mean).dot(unit_v)))

	results = []
	for i, shape in enumerate(shapes):
		result = []
		for j,measurement in enumerate(measures):
			result.append(measurement.execute(shape))
		results.append(result)


	
	writeReport('reportTest.xlsx', measures, results)
