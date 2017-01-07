#! /usr/bin/python

USE_CACHE = False

import os
import sys
import timeit
import numpy as np
import time
import argparse

import config.defaults
import src.spherical_harmonics as sh
from src.vertexarea import getVertexAreas
from src.shapes import Sphere
from src.reporting import writeCSV

from cycler import cycler
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as mplcm


def processSphere(filePath, nominalRadius, Lmax, verbose):
	if verbose:
		print 'Processing ' + os.path.split(filePath)[-1][:-4]

	cacheFileName = filePath[:-4] + '_cached_' + str(Lmax)
	if os.path.isfile(cacheFileName + '.npy') and USE_CACHE:
		finalYs = np.load(cacheFileName + '.npy')
	else:
		sphere = Sphere(filePath, nominalRadius, False)

		sphericalCoordinates = sh.getSphericalCoordinates(sphere.vertices, sphere.centerPoint)

		if verbose:
			print 'Calculating areas'
		vertexAreas = getVertexAreas(sphere.faces, sphere.vertices)
		if verbose:
			print 'Finished calculating areas'

		finalYs, coefficients = sh.simple_transform(sphericalCoordinates, Lmax, vertexAreas)
		if verbose:
			print 'Coefficients:'
			for l in range(Lmax):
				print l
				print coefficients[l**2:(l**2 + 2*l + 1)]
		np.save(cacheFileName, finalYs)

	return finalYs

def checkOBJpath(path):
	if not os.path.isdir(path) and (not os.path.isfile(path) and path.endswith('.obj')):
		raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(directory))
	return path

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Generates a spreadsheet')
	parser.add_argument("path", help="The folder or the obj-file", type=checkOBJpath)
	parser.add_argument("frequencies", help="The number of frequencies to analyze", type=int)
	parser.add_argument("--radius", help="The nominal radius", type=float)
	parser.add_argument("--verbose", help="Show debug information", action='store_true')
	parser.add_argument("--csv", help="Instead of plot, write csv file", action='store_true')
	parser.add_argument("-o", "--output", help="Path for csv file", default='spherical_harmonics.csv', type=str)
	args = parser.parse_args()

	if args.radius is None:
		args.radius = config.defaults.nominalRadius
	
	if args.verbose:
		print 'Using radius ' + str(args.radius)

	if args.path.endswith('.obj'):
		folderPath = ''
		files = [args.path]
	else:
		folderPath = args.path
		files = [fileName for fileName in os.listdir(folderPath) if fileName.endswith('.obj')]

	yss = []
	if not args.csv:
		numColors = 0
		for fileName in files:
			if fileName.endswith('.obj'):
				numColors += 1

		cm = plt.get_cmap('gist_rainbow')

		fig = plt.figure()
		ax = fig.add_subplot(111)
		cNorm  = mcolors.Normalize(vmin=0, vmax=numColors-1)
		scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)

		colors = [scalarMap.to_rgba(i) for i in range(numColors)]
		lineStyles = ['-', '--', ':', '-.']
		lineWidths = [2,3]

		paddedLineStyles = np.tile(lineStyles, len(colors)/len(lineStyles))
		paddedLineWidths = np.tile(lineWidths, len(colors)/len(lineWidths))
		if len(paddedLineStyles) < len(colors):
			paddedLineStyles = np.hstack((paddedLineStyles,lineStyles[0:len(colors)-len(paddedLineStyles)]))
		if len(paddedLineWidths) < len(colors):
			paddedLineWidths = np.hstack((paddedLineWidths,lineWidths[0:len(colors)-len(paddedLineWidths)]))

		ax.set_prop_cycle(cycler('color', colors) + cycler('linestyle', paddedLineStyles) + cycler('linewidth', paddedLineWidths))

	for fileName in files:
		start_time = timeit.default_timer()
		ys = processSphere(os.path.join(folderPath,fileName), args.radius, args.frequencies, args.verbose)
		if args.verbose:
			print 'Total: ' + str(np.sum(ys))
			print('Time taken: ' + str(timeit.default_timer() - start_time))

		yss.append(ys.tolist())
	if args.csv:
		if len(yss) == 1:
			yss = yss[0]
			yss = [[val] for val in yss]
		else:
			frequencies = np.arange(args.frequencies+1)
			print frequencies.shape
			yss = np.array(yss).T
			yss = np.insert(yss,0,frequencies, axis=1)

		writeCSV(args.output, yss)
	else:
		for ys,fileName in zip(yss,files):
			xa = np.arange(0, len(ys))
			plt.plot(xa, ys, label=fileName[0:-4])

		plt.legend(fontsize=9)
		plt.grid(True)
		plt.show()
