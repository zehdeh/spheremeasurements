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

if __name__ == '__main__':
	
	if len(sys.argv) < 2:
		print 'Please provide as an argument the directory where the OBJ-files are located and the radius to be used'
		sys.exit(1)

	shapes = []

	for fileName in os.listdir(sys.argv[1]):
		if fileName.endswith('.obj'):
			filePath = sys.argv[1] + '/' + fileName
			print 'Loading mesh ' + fileName
			shapes.append(Sphere(filePath, *sys.argv[2:]))

	
	measures = getMeasures(Sphere)

	centerPoints = [s.centerPoint for s in shapes]
	mean, unit_v, projectedPoints, pointDistances = scanLineFit(centerPoints)

	measures.append(Measure('Distance from cam', lambda x: (x.centerPoint - mean).dot(unit_v)))

	results = []
	for i, shape in enumerate(shapes):
		result = []
		for j,measurement in enumerate(measures):
			result.append(measurement.execute(shape))
		results.append(result)


	
	writeReport('reportTest.xlsx', measures, results)
