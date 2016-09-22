#! /usr/bin/python

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

def centerModel(vertices):
	avgs = [np.mean(x) for x in vertices.T]
	vertices -= avgs
	return avgs

def loadOBJ(fileName):
	"""Loads a Wavefront OBJ file. """
	vertices = []
	normals = []
	texcoords = []
	faces = []

	material = None
	for line in open(fileName, "r"):
		if line.startswith('#'): continue
		values = line.split()
		if not values: continue
		if values[0] == 'v':
			 v = map(float, values[1:4])
			 vertices.append(v)
		elif values[0] == 'vn':
			 v = map(float, values[1:4])
			 normals.append(v)
		elif values[0] == 'vt':
			 texcoords.append(map(float, values[1:3]))
		elif values[0] == 'f':
			 face = []
			 texcoords = []
			 norms = []
			 for v in values[1:]:
				if '//' in v:
					glue = '//'
				else:
					glue = '/'
				w = v.split(glue)
				face.append(int(w[0]) - 1)
				if len(w) >= 2 and len(w[1]) > 0:
					texcoords.append(int(w[1]))
				else:
					texcoords.append(0)
					if len(w) >= 3 and len(w[2]) > 0:
						norms.append(int(w[2]))
					else:
						norms.append(0)
			 faces.append(face)
	return np.asarray(vertices), np.asarray(faces), np.asarray(normals)

def writeOBJ(fileName, vertices, faces, normals):
	with open(fileName, 'w') as f:
		f.write("# OBJ file\n")
		for v in vertices:
			f.write('v ' +' '.join([format(x,'.4f') for x in v]))
			f.write("\n")
		for p in faces:
			f.write("f")
			for i in p:
				f.write(" %d" % (i + 1))
			f.write("\n")

def removeVerticesByCondition(condition, vertices, faces):
	mask = condition(vertices)
	indices = [i for i,x in enumerate(zip(vertices,mask)) if x[1]]
	indicesExclude = np.asarray([i for i,j in enumerate(vertices) if i not in indices])

	vertices = vertices[indices]
	faces = [f for f in faces if not bool(set(indicesExclude) & set(f))]

	for i,f in enumerate(faces):
		faces[i] = [x - len(indicesExclude[x > indicesExclude]) for x in f]

	return vertices, np.asarray(faces)

def removeIsolatedVertices(vertices, faces):
	referencedVertices = [i for f in faces for i in f]
	indices = set([i for i,v in enumerate(vertices)])
	isolatedIndices = np.asarray(list(indices - set(referencedVertices)))
	vertices = vertices[list(set(referencedVertices))]

	for i,f in enumerate(faces):
		faces[i] = [x - len(isolatedIndices[x > isolatedIndices]) for x in f]
	
	return vertices, faces

def houghTransformation(vertices, faces):
	verticesZ = vertices
	verticesZ.T[2] = 0
	print verticesZ[0]

	minX = np.min(verticesZ.T[0])
	maxX = np.max(verticesZ.T[0])
	sizeX = maxX - minX
	minY = np.min(verticesZ.T[1])
	maxY = np.max(verticesZ.T[1])
	sizeY = maxY - minY

	grid = np.zeros((sizeX,sizeY))
	plt.scatter(verticesZ.T[0], vertices.T[1])
	plt.show()

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

			#uv = np.asarray([0,1,0])
			#uv = uv.T / np.linalg.norm(uv)
			#plane = uv.tolist() + [-200]
			#condition = lambda x: pointPlaneDistance(plane, x.T) < 0
			#centerPoint = np.asarray([270,1460,-200])
			#centerPoint = np.asarray([0,0,0])
			#condition = lambda x: np.linalg.norm(centerPoint - x, axis=1) > 1340
			#vertices, faces = removeVerticesByCondition(condition, vertices, faces)
			#vertices, faces = removeIsolatedVertices(vertices, faces)
			#offset = centerModel(vertices)
			houghTransformation(vertices, faces)

			if sys.argv[2].endswith('.obj'):
				writeOBJ(sys.argv[2], vertices, faces, normals)
			else:
				outputFile = sys.argv[2] + '/' + fileName
				writeOBJ(outputFile, vertices, faces, normals)


