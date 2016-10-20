#! /usr/bin/python


import os
import sys
import timeit
import numpy as np
import time
import pymesh

import Queue
import threading

from math import isnan, sin
from scipy.special import sph_harm
from src.OBJIO import loadOBJ
from src.fitting import fitSphere
from src.mesh import Mesh
from scipy.optimize import leastsq

from cycler import cycler
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def getSphericalCoordinates(vertices, centerPoint):
	vertices = vertices - centerPoint
	r = np.sqrt(vertices.T[0]**2 + vertices.T[1]**2 + vertices.T[2]**2)
	theta = np.arccos(vertices.T[2] / r)
	phi = np.arctan(vertices.T[1]/vertices.T[0])
	return np.array([phi, theta, r])

def getCartesianCoordinates(theta, phi, r, centerPoint):
	x = r*np.sin(theta)*np.cos(phi)
	y = r*np.sin(theta)*np.sin(phi)
	z = r*np.cos(theta)

	coords = np.array([x,y,z])
	coords = coords.T + centerPoint

	return coords

def residuals(aa, Yval, rs):
	return (Yval.dot(aa) - rs)

def getAMatrix(theta, phi, l):
	n = len(phi)

	A = np.zeros((n, 2*l+1))
	for i in range(n):
		A[i] = [sph_harm(l, j-l, theta[i], phi[i]).real for j in range(0, 2*l+1)]
	A = np.nan_to_num(A)
	
	return A
def sirf(sphericalCoordinates, Lmax):
	phi, theta, r = sphericalCoordinates

	s = 5
	A = getAMatrix(theta, phi, 0)
	for i in range(1,s):
		A = np.concatenate((A,getAMatrix(theta, phi, i)), axis=1)
	
	b = np.zeros((Lmax+1)**2)
	b[0:A.shape[1]], flag = leastsq(residuals, b[0:A.shape[1]], args=(A,r))

	res = r - A.dot(b[0:A.shape[1]])
	for l in range(s, Lmax+1):
		Al = getAMatrix(theta, phi, l)
		lowerBound = l**2
		upperBound = (l+1)**2

		b[lowerBound:upperBound], flag = leastsq(residuals, b[lowerBound:upperBound], args=(Al,res))
		res = res - Al.dot(b[lowerBound:upperBound])
		A = np.concatenate((A,Al), axis=1)


	ys = np.zeros(Lmax+1)
	for i in range(Lmax+1):
		lowerBound = i**2
		upperBound = (i+1)**2
		ys[i] = np.linalg.norm(b[lowerBound:upperBound], ord=2)
	
	return b, ys, A

def nlsf(sphericalCoordinates, LMax):
	phi, theta, r = sphericalCoordinates

	n = len(phi)
	k = (Lmax + 1)**2
	Y = np.zeros((n,k))
	
	kToLM = np.zeros((k,2))
	for i in range(n):
		for l in range(Lmax+1):
			for m in range(-l,l+1):
				j = l**2 + l + m + 1
				kToLM[j-1,:] = [l,m]
				s = sph_harm(l, m, theta[i], phi[i]).real
				#print 'l ' + str(l)
				#print 'm ' + str(m)
				if not isnan(s):
					Y[i,j-1] = s

	a0 = np.zeros(k)
	a1, flag = leastsq(residuals, a0, args=(Y,r))

	ys = np.zeros(Lmax+1)
	LMMatrix = np.zeros((k,k))
	for j,a in enumerate(a1):
		l,m = kToLM[j]
		LMMatrix[int(l),int(m+l)] = a
		#ys[int(l)]  = np.linalg.norm(LMMatrix[int(l)], ord=2)
		ys[int(l)] = np.sum(LMMatrix[int(l)])

	r_approx = Y.dot(a1)

	return a1, ys, r_approx

def loadOBJwithSphericalCoordinates(fileName):
	vertices, faces, normals = loadOBJ(fileName)
	bounds = Mesh(vertices.T, faces, normals).getBounds()
	p0 = [bounds[0][0],bounds[1][0],bounds[2][0],150]
	centerPoint, radius = fitSphere(vertices, p0, 150, bounds)

	mesh = pymesh.form_mesh(vertices, faces)
	mesh.add_attribute('vertex_area')
	vertexAreas = mesh.get_attribute('vertex_area')

	sphericalCoordinates = getSphericalCoordinates(vertices, centerPoint)
	#phi, theta, r = sphericalCoordinates
	#sphericalCoordinates = np.array([phi, theta, r-radius])

	return sphericalCoordinates, vertexAreas

def processSphere(filePath):
	Lmax = int(sys.argv[2])
	sphericalCoordinates, vertexAreas = loadOBJwithSphericalCoordinates(filePath)

	phi, theta, r = sphericalCoordinates
	vertexAreas = vertexAreas / np.max(vertexAreas)
	#vertexAreas = 1/vertexAreas

	r = r*(vertexAreas)
	sphericalCoordinates = [phi, theta, r]
	noPasses = 2

	finalYs = np.zeros(Lmax + 1)
	for i in range(noPasses):
		print "Pass " + str(i)
		a, ys, A = sirf(sphericalCoordinates, Lmax)
		finalYs += ys
		phi, theta, r = sphericalCoordinates

		r = r - A.dot(a)
		rmse = np.linalg.norm(r) / np.sqrt(len(r))
		print "RMSE: " + str(rmse)
		sphericalCoordinates = np.array([phi, theta, r])

		
	return finalYs

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print "Please specify OBJPATH NOFREQUENCIES"
		sys.exit(0)


	if sys.argv[1].endswith('.obj'):
		folderPath = ''
		files = [sys.argv[1]]
	else:
		folderPath = sys.argv[1]
		files = os.listdir(folderPath)
	numColors = len(files)

	cm = plt.get_cmap('gist_rainbow')

	fig = plt.figure()
	ax = fig.add_subplot(111)
	#colors = [cm(1.*i/numColors) for i in range(numColors)]
	#lineStyle = np.tile(['-', '--', ':', '-.'], len(colors)/4)
	#lineWidth = np.tile([2,3], len(colors)/2)

	#ax.set_prop_cycle(cycler('color', colors) + cycler('linestyle', lineStyle) + cycler('linewidth', lineWidth))
	#ax.set_prop_cycle(cycler('color', colors) + cycler('linewidth', lineWidth))
	'''
	threads = list()
	for i in range(10):
		thread = spThread('thread-' + str(i), workQueue, queueLock)
		thread.start()
		threads.append(thread)
		'''
	for fileName in files:
		if fileName.endswith('.obj'):
			#workQueue.put(folderPath + fileName)
			#sphericalCoordinates = loadOBJwithSphericalCoordinates(folderPath + fileName)

			start_time = timeit.default_timer()
			ys = processSphere(folderPath + fileName)
			print(timeit.default_timer() - start_time)

			xa = np.arange(0, len(ys))
			plt.plot(xa, ys, label=fileName[0:-4])

			#rmse = np.linalg.norm(r - r_approx) / np.sqrt(len(a1))

		plt.legend()
	plt.grid(True)
	plt.gca().set_yscale('log')
	plt.show()
