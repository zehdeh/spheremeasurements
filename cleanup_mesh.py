#! /usr/bin/python

import sys
import numpy as np

def centerModel(vertices):
	avgs = [np.mean(x) for x in vertices.T]
	vertices -= avgs

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
				  w = v.split('/')
				  face.append(int(w[0]))
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
				f.write(" %d" % i)
			f.write("\n")

def removeVerticesByCondition(condition, vertices, faces):
	mask = vertices[condition(vertices)]
	indices = np.where(mask)
	indicesExclude = np.where([not m for m in mask])
	print indices
	print indicesExclude


	vertices = vertices[indices]
	facesIncluded = []
	for f in faces:
		includeF = True
		for i in indicesExclude:
			if i in f:
				includeF = False
		if includeF:
			facesIncluded.append(f)
	print len(faces)
	print len(facesIncluded)
	return vertices, np.asarray(facesIncluded)

def pointPlaneDistance(point, vertices):
	plane_xyz = point[0:3]
	distance = (plane_xyz*vertices.T).sum(axis=1) + point[3]
	return distance / np.linalg.norm(plane_xyz)

if __name__ == '__main__':
	fileName = sys.argv[1]
	vertices, faces, normals = loadOBJ(fileName)

	centerModel(vertices)

	plane = [-1, 0, 0, 0]
	condition = lambda x: pointPlaneDistance(plane, x.T) > 10000
	vertices, faces = removeVerticesByCondition(condition, vertices, faces)

	outputFile = sys.argv[2]
	writeOBJ(outputFile, vertices, faces, normals)


