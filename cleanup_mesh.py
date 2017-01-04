#! /usr/bin/python

import sys
import os
import numpy as np
from math import floor, ceil
import matplotlib.pyplot as plt
from src.OBJIO import loadOBJ, writeOBJ, getVTKMesh
from opendr.topology import get_vert_connectivity
from scipy.sparse.csgraph import connected_components
from src.fitting import calculateMeanCurvature
import scipy
import networkx as nx
from src.vertexarea import getFaceAngles, getFaceArea

def centerModel(vertices):
	avgs = [np.mean(x) for x in vertices.T]
	vertices -= avgs
	return avgs

def removeSmallIsolatedComponents(vertices, faces, normals):
	cnct = get_vert_connectivity(vertices, faces)
	nbrs = [np.nonzero(np.array(cnct[:,i].todense()))[0] for i in range(cnct.shape[1])]
	A = cnct.toarray()
	A[A > 0] = 1

	G = nx.from_numpy_matrix(A)
	connectedComponents = sorted(nx.connected_components(G), key=len, reverse=True)
	if len(connectedComponents) > 1:
		indicesToDelete = np.array([j for c in connectedComponents[1:] for j in c])
		vertices = np.delete(vertices, indicesToDelete, axis=0)
		normals = np.delete(normals, indicesToDelete, axis=0)

		faceHasBadVs = np.all(np.in1d(faces.ravel(), indicesToDelete, invert=True).reshape(faces.shape), axis=1)
		faces = faces[np.where(faceHasBadVs)]
		faces = np.array([f - len(indicesToDelete[f > indicesToDelete]) for f in faces.ravel()]).reshape(faces.shape)

	return vertices, faces, normals

def removeVerticesByCondition(condition, vertices, faces, normals):
	mask = condition(vertices)
	indices = [i for i,x in enumerate(zip(vertices,mask)) if x[1]]

	return removeVerticesByIndex(vertices, faces, normals, indices)


def removeIsolatedVertices(vertices, faces, normals):
	referencedVertices = [i for f in faces for i in f]
	indices = set([i for i,v in enumerate(vertices)])
	isolatedIndices = np.asarray(list(indices - set(referencedVertices)))
	vertices = vertices[list(set(referencedVertices))]
	normals = normals[list(set(referencedVertices))]

	faces = np.array([f - len(isolatedIndices[f > isolatedIndices]) for f in faces.ravel()]).reshape(faces.shape)
	
	return vertices, faces, normals


def removeVerticesByIndex(vertices, faces, normals, indices):
	indicesExclude = np.asarray([i for i,j in enumerate(vertices) if i not in indices])

	vertices = vertices[indices]
	normals = normals[indices]
	faceHasBadVs = np.all(np.in1d(faces.ravel(), indicesExclude, invert=True).reshape(faces.shape), axis=1)
	faces = faces[np.where(faceHasBadVs)]

	faces = np.array([f - len(indicesExclude[f > indicesExclude]) for f in faces.ravel()]).reshape(faces.shape)

	return vertices, np.asarray(faces), normals

def removePointsWithExtremeCurvature(vertices, faces, normals):
	curvatureThreshold = 3
	polyData = getVTKMesh(vertices,faces,normals)
	curvature = calculateMeanCurvature(polyData)
	indices = np.where(curvature < curvatureThreshold)[0]
	indicesExclude = np.asarray([i for i,j in enumerate(vertices) if i not in indices])

	vertices = vertices[indices]
	normals = normals[indices]
	faceHasBadVs = np.all(np.in1d(faces.ravel(), indicesExclude, invert=True).reshape(faces.shape), axis=1)
	faces = faces[np.where(faceHasBadVs)]

	faces = np.array([f - len(indicesExclude[f > indicesExclude]) for f in faces.ravel()]).reshape(faces.shape)

	return vertices, np.asarray(faces), normals

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

def processMesh(fileName, folderPath):
		print 'Processing ' + fileName + '...'
		vertices, faces, normals = loadOBJ(folderPath + '/' + fileName)

		radiusNominal = 80.065605

		sphereCenter = houghTransformation(vertices, faces, normals, radiusNominal)
		condition = lambda x: np.linalg.norm(sphereCenter - x, axis=1) < (radiusNominal*1.025)
		vertices, faces, normals = removeVerticesByCondition(condition, vertices, faces, normals)

		#vertices, faces, normals = removeSmallIsolatedComponents(vertices, faces, normals)

		'''
		faceLargestAngles = np.zeros(faces.shape[0])
		for i,f in enumerate(faces):
			a,b,c = getFaceAngles(f, vertices)
			faceLargestAngles[i] = max(a,b,c)
		
		print faces[np.where(faceLargestAngles>np.pi - 0.2)]
		faces = faces[np.where(faceLargestAngles<np.pi - 0.2)]
		'''
		vertices, faces, normals = removeIsolatedVertices(vertices, faces, normals)

		vertices, faces, normals = removePointsWithExtremeCurvature(vertices, faces, normals)
		
		if sys.argv[2].endswith('.obj'):
			writeOBJ(sys.argv[2], vertices, faces, normals)
		else:
			outputFile = sys.argv[2] + '/' + fileName
			writeOBJ(outputFile, vertices, faces, normals)

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
			processMesh(fileName, folderPath)
	


