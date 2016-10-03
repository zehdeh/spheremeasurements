#! /usr/bin/python

import sys
import numpy as np
from math import isnan
from scipy.special import sph_harm
from src.OBJIO import loadOBJ
from src.fitting import fitSphere
from src.mesh import Mesh
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

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

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print "Please specify OBJPATH NOFREQUENCIES"
		sys.exit(0)
	vertices, faces, normals = loadOBJ(sys.argv[1])

	bounds = Mesh(vertices.T, faces, normals).getBounds()
	p0 = [bounds[0][0],bounds[1][0],bounds[2][0],150]

	centerPoint, radius = fitSphere(vertices, p0, 150, bounds)
	
	#fig = plt.figure()
	#ax = fig.add_subplot(111, projection='3d')
	#ax.scatter(vertices.T[0], vertices.T[1], vertices.T[2])
	#ax.scatter(centerPoint[0], centerPoint[1], centerPoint[2], color='r')
	#plt.show()

	# Create a sphere
	pi = np.pi
	cos = np.cos
	sin = np.sin
	#phi = np.linspace(0,pi,101)
	#theta = np.linspace(-pi,pi,101)
	phi, theta, r = getSphericalCoordinates(vertices, centerPoint)

	x = r * sin(phi) * cos(theta)
	y = r * sin(phi) * sin(theta)
	z = r * cos(phi)

	Lmax = int(sys.argv[2])
	n = len(phi)

	k = (Lmax + 1)**2
	Y = np.zeros((n,k))
	# Represent spherical harmonics on the surface of the sphere
	for i in range(n):
		for l in range(Lmax+1):
			for m in range(-l,l+1):
				j = l**2 + l + m + 1
				s = sph_harm(l, m, theta[i], phi[i]).real
				#print 'l ' + str(l)
				#print 'm ' + str(m)
				if not isnan(s):
					Y[i,j-1] = s

	a0 = np.zeros(k)
	a1, flag = leastsq(residuals, a0, args=(Y,r))

	r_approx = Y.dot(a1)

	vertices2 = getCartesianCoordinates(theta, phi, r_approx, centerPoint)
	rmse = np.linalg.norm(r - r_approx) / np.sqrt(len(a1))

	xplot = np.arange(0, len(a1), 1)
	fig = plt.figure()
	#plt.plot(xplot, a1)
	a1 = np.fabs(a1)
	bins = np.arange(a1.min(), a1.max()-1,1)
	#bins = np.logspace(a1.min(), a1.max()-1,1000)
	#plt.xticks(bins, ["2^%s" % i for i in bins])
	plt.hist(a1, bins=bins)
	#plt.hist(a1, bins=np.arange(a1.min(), a1.max()-1))
	#plt.gca().set_xscale('log')
	#plt.gca().set_yscale('log')
	plt.show()


	#x2plot = np.arange(0, len(rmses), 1)
	#plt.plot(x2plot, rmses)
	plt.show()
