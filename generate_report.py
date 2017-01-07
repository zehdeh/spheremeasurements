#! /usr/bin/python

import os
import sys
import time
import importlib
import opendr
import numpy as np
import argparse
import os

import config.defaults
from src.reporting import writeReport
from src.measure import Measure, getMeasures
from src.shapes import Sphere
from src.utils import checkDir

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Generates a spreadsheet')
	parser.add_argument("folder", help="The folder with the OBJ files", type=checkDir)
	parser.add_argument("--radius", help="The nominal radius to be used", type=float)
	parser.add_argument("--fit-radius", help="Fit the radius instead of using the nominal value", action='store_true', default=False)
	parser.add_argument("-o", "--output", help="The name of the file that is to be written", type=str, default='report')
	parser.add_argument("--csv", help="Generate a raw CSV file instead of an excel sheet", action='store_true', default=False)
	parser.add_argument("--calibration", help="If you are interested in camera distance, provide a calibration folder", type=checkDir)
	parser.add_argument("--camera", help="The camera to be used for distance measurement", type=int)
	parser.add_argument("--extract-num", help="Extract the number included in the filename", action='store_true')
	parser.add_argument("--verbose", help="Show debug information", action='store_true')

	args = parser.parse_args()

	if args.radius is None:
		args.radius = config.defaults.nominalRadius
	
	if args.verbose:
		print 'Using radius ' + str(args.radius)
	
	spheres = []
	for fileName in os.listdir(args.folder):
		if fileName.endswith('.obj'):
			filePath = sys.argv[1] + '/' + fileName
			if args.verbose:
				print 'Loading file ' + fileName
			spheres.append(Sphere(filePath, args.radius, args.fit_radius))
	
	if len(spheres) == 0:
		print 'No meshes found!'
		sys.exit(0)

	measures = getMeasures(args.calibration, args.camera, args.extract_num)

	sortColumn = -1
	results = np.zeros((len(spheres),len(measures)))
	for i, sphere in enumerate(spheres):
		for j,measurement in enumerate(measures):
			results[i,j] = measurement.execute(sphere)
			if measurement.sort:
				sortColumn = j

	if sortColumn != -1:
		sortedIndices = np.argsort(results.T[sortColumn])
		results = results[sortedIndices]

		spheres = [spheres[i] for i in sortedIndices]
	results = results.tolist()
		
	for i, sphere in enumerate(spheres):
		results[i] = [sphere.fileName] + results[i]
	
	if args.output == 'report':
		outputName = args.output + '.xlsx' if not args.csv else args.output + '.csv'
	else:
		outputName = args.output

	if args.verbose:
		print 'Writing to ' + os.path.abspath(outputName)
	writeReport(outputName, measures, results, args.csv)
