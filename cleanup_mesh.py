#! /usr/bin/python

import sys
import os
import numpy as np
from math import floor, ceil
import matplotlib.pyplot as plt
from src.OBJIO import loadOBJ, writeOBJ

def centerModel(vertices):
	avgs = [np.mean(x) for x in vertices.T]
	vertices -= avgs
	return avgs


def removeVerticesByCondition(condition, vertices, faces, normals):
	mask = condition(vertices)
	indices = [i for i,x in enumerate(zip(vertices,mask)) if x[1]]
	indicesExclude = np.asarray([i for i,j in enumerate(vertices) if i not in indices])

	vertices = vertices[indices]
	normals = normals[indices]
	faces = [f for f in faces if not bool(set(indicesExclude) & set(f))]

	for i,f in enumerate(faces):
		faces[i] = [x - len(indicesExclude[x > indicesExclude]) for x in f]

	return vertices, np.asarray(faces), normals

def removeIsolatedVertices(vertices, faces, normals):
	referencedVertices = [i for f in faces for i in f]
	indices = set([i for i,v in enumerate(vertices)])
	isolatedIndices = np.asarray(list(indices - set(referencedVertices)))
	vertices = vertices[list(set(referencedVertices))]
	normals = normals[list(set(referencedVertices))]

	for i,f in enumerate(faces):
		faces[i] = [x - len(isolatedIndices[x > isolatedIndices]) for x in f]
	
	return vertices, faces, normals

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

def pointPlaneDistance(point, vertices):
	plane_xyz = point[0:3]
	distance = (plane_xyz*vertices.T).sum(axis=1) + point[3]
	return distance / np.linalg.norm(plane_xyz)

if __name__ == '__main__':
	if sys.argv[1].endswith('.obj'):
		path = os.path.split(sys.argv[1])
		files = [path[1]]
		folderPath = path[0]
		if folderPath == '':
			folderPath = './'
	else:
		folderPath = sys.argv[1]
		files = os.listdir(folderPath)
	
	for fileName in files:
		if fileName.endswith('.obj'):
			print 'Processing ' + fileName + '...'
			vertices, faces, normals = loadOBJ(folderPath + '/' + fileName)

			vertices, faces, normals = removeIsolatedVertices(vertices, faces, normals)
			#uv = np.asarray([0,1,0])
			#uv = uv.T / np.linalg.norm(uv)
			#plane = uv.tolist() + [-200]
			#condition = lambda x: pointPlaneDistance(plane, x.T) < 0
			#centerPoint = np.asarray([270,1460,-200])
			#centerPoint = np.asarray([0,0,0])
			#condition = lambda x: np.linalg.norm(centerPoint - x, axis=1) > 1340
			#vertices, faces = removeVerticesByCondition(condition, vertices, faces)
			#offset = centerModel(vertices)
			sphereCenter = houghTransformation(vertices, faces, normals, 150)
			condition = lambda x: np.linalg.norm(sphereCenter - x, axis=1) < 160
			vertices, faces, normals = removeVerticesByCondition(condition, vertices, faces, normals)

			if sys.argv[2].endswith('.obj'):
				writeOBJ(sys.argv[2], vertices, faces, normals)
			else:
				outputFile = sys.argv[2] + '/' + fileName
				writeOBJ(outputFile, vertices, faces, normals)


