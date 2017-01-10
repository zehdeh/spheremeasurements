#! /usr/bin/python

import sys
import os
import argparse
import numpy as np
from math import floor,ceil

import config.defaults
from src.OBJIO import loadOBJ, writeOBJ, getVTKMesh
from src.utils import checkFile

def houghTransformation(vertices, faces, normals, radius):
	winners = np.zeros((2,3))
	for i,dim in enumerate([1,2]):
		otherDims = [0,1,2]
		otherDims.remove(dim)

		lowerBounds = [np.min(v) for v in vertices.T[otherDims]]
		upperBounds = [np.max(v) for v in vertices.T[otherDims]]
		sizes = [upperBounds[0] - lowerBounds[0], upperBounds[1] - lowerBounds[1]]

		minX = np.min(vertices.T[0])
		maxX = np.max(vertices.T[0])
		minY = np.min(vertices.T[1])
		maxY = np.max(vertices.T[1])

		centers = []
		for v,n in zip(vertices, normals):
			if np.count_nonzero(n) > 0:
				n = n/np.linalg.norm(n)
			n[dim] = 0
			n = -n
			c = v+n*radius
			if c[otherDims[0]] > lowerBounds[0] and c[otherDims[0]] < upperBounds[0] and c[otherDims[1]] > lowerBounds[1] and c[otherDims[1]] < upperBounds[1]:
				centers.append(c)
		centers = np.asarray(centers)

		grid = np.zeros((int(ceil(sizes[0])),int(ceil(sizes[1]))))
		for c in centers:
			a = c[otherDims[0]] - lowerBounds[0]
			b = c[otherDims[1]] - lowerBounds[1]
			j = int(floor(a))
			k = int(floor(b))
			grid[j,k] += 1
		
		winnerIndex = np.unravel_index(np.argmax(grid), grid.shape)
		winners[i,otherDims[0]] = lowerBounds[0] + winnerIndex[0]
		winners[i,otherDims[1]] = lowerBounds[1] + winnerIndex[1]

		#plt.scatter(vertices.T[otherDims[0]], vertices.T[otherDims[1]])
		#plt.scatter(centers.T[0], centers.T[1], color='r')
		#plt.scatter(winners[i,otherDims[0]], winners[i,otherDims[1]], color='r')
		#plt.show()
	return np.array([np.mean(winners.T[0]), winners[1,1], winners[0,2]])

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Generates a spreadsheet')
	parser.add_argument("file", help="The OBJ file to clean", type=lambda x: checkFile(x,'.obj'))
	parser.add_argument("--verbose", help="Show debug information", action='store_true')
	args = parser.parse_args()

	nominalRadius = config.defaults.nominalRadius

	vertices, faces, normals = loadOBJ(args.file)
	print houghTransformation(vertices, faces, normals, nominalRadius)

	print vertices.shape
