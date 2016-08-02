#! /usr/bin/python

import sys
import math, random
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def generateSphere(samples=1,radius=1,randomize=True):
	rnd = 1.
	if randomize:
		rnd = random.random() * samples
	
	points = []
	offset = 2./samples
	increment = math.pi * (3. - math.sqrt(5.))

	for i in range(samples):
		y = ((i * offset) - 1) + (offset / 2)
		r = math.sqrt(1 - pow(y,2))

		phi = ((i + rnd) % samples) * increment

		x = math.cos(phi) * r * radius
		z = math.sin(phi) * r * radius
		points.append([x,y,z])
	
	return np.asmatrix(points)


def loadOBJ(filename):
	numVerts = 0
	vertices = np.array([])
	normals = np.array([])
	verts = []
	norms = []
	vertsOut = []
	normsOut = []
	for line in open(filename, "r"):
		vals = line.split()
		if len(vals) == 0:
			continue
		if vals[0] == "v":
			v = map(float, vals[1:4])
			verts.append(v)
		if vals[0] == "vn":
			n = map(float, vals[1:4])
			norms.append(n)
		'''
		if vals[0] == "f":
			for f in vals[1:]:
				w = f.split("/")
				# OBJ Files are 1-indexed so we must subtract 1 below
				vertsOut.append(list(verts[int(w[0])-1]))
				normsOut.append(list(norms[int(w[2])-1]))
				numVerts += 1
		'''
	return np.asmatrix(verts), normsOut

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print "Please specify a file to open"
		sys.exit(1)

	vertices,norms = loadOBJ(sys.argv[1])
	x,y,z = vertices.T
	x = x[0].tolist()[0]
	y = y[0].tolist()[0]
	z = z[0].tolist()[0]
	
	fig = plt.figure()
	ax = fig.add_subplot(111,projection='3d')
	ax.scatter(x,y,z)
	plt.show()
