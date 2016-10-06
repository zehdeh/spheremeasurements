#! /usr/bin/python



import sys
import timeit
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

def getAMatrix(theta, phi, l):
	n = len(phi)
	A = np.zeros((n,2*l+1))
	for i in range(n):
		for j in range(0,2*l+1):
			s = sph_harm(l, j-l, theta[i], phi[i]).real
			if not isnan(s):
				A[i, j] = s
	
	return A
def sirf(sphericalCoordinates, LMax):
	phi, theta, r = sphericalCoordinates

	s = 5
	A = getAMatrix(theta, phi, 0)
	for i in range(1,s):
		A = np.concatenate((A,getAMatrix(theta, phi, i)), axis=1)
	
	b0 = np.zeros(A.shape[1])
	b, flag = leastsq(residuals, b0, args=(A,r))

	res = r - A.dot(b)
	for l in range(s, Lmax+1):
		Al = getAMatrix(theta, phi, l)
		bl0 = np.zeros(Al.shape[1])

		bl, flag = leastsq(residuals, bl0, args=(Al,res))
		res = res - Al.dot(bl)

		lowerBound = l**2
		upperBound = (l+1)**2
		print lowerBound
		print upperBound
		print len(bl)
		b = np.concatenate((b, bl))

	LMatrix = np.zeros((Lmax+1, Lmax*2+1))
	ys = np.zeros(Lmax+1)
	for i in range(Lmax+1):
		lowerBound = i**2
		upperBound = (i+1)**2
		LMatrix[i, 0:(upperBound-lowerBound)] = b[lowerBound:upperBound]
		ys[i] = np.linalg.norm(LMatrix[i, :], ord=2)
		#print 'from ' + str((i)**2)
		#print 'to ' + str((i+1)**2)
	
	return b, ys

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
		ys[int(l)]  = np.linalg.norm(LMMatrix[int(l)], ord=2)

	r_approx = Y.dot(a1)

	return a1, ys, r_approx

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
	sphericalCoordinates = getSphericalCoordinates(vertices, centerPoint)
	Lmax = int(sys.argv[2])
	phi, theta, r = sphericalCoordinates

	#a1, ys, r_approx = nlsf(sphericalCoordinates, Lmax)
	
	#a2, ys2 = sirf(sphericalCoordinates, Lmax)
	start_time = timeit.default_timer()
	nlsf(sphericalCoordinates, Lmax)
	print(timeit.default_timer() - start_time)

	start_time = timeit.default_timer()
	sirf(sphericalCoordinates, Lmax)
	print(timeit.default_timer() - start_time)

	sys.exit(0)

	#vertices2 = getCartesianCoordinates(theta, phi, r_approx, centerPoint)
	#rmse = np.linalg.norm(r - r_approx) / np.sqrt(len(a1))

	xplot = np.arange(0, len(a2), 1)
	fig = plt.figure()
	#plt.plot(xplot, a1)
	#a1 = np.fabs(a1)
	#a1 = np.sort(a1)[::-1]
	#a1 = a1[a1 > 0]
	#print a1.max()

	xa = np.arange(0, len(ys2))
	#plt.plot(xa, ys)
	plt.plot(xa, ys2, color='r')
	plt.grid(True)
	plt.title(sys.argv[1])
	plt.gca().set_yscale('log')
	plt.show()

	#bins = np.arange(a1.min(),a1.max(),len(a1))
	#bins = 10**np.linspace(np.log10(0.00001), np.log10(600),50)
	#bins = np.logspace(0.00000001, 0.000001, 50)
	#plt.xticks(bins, ["2^%s" % i for i in bins])
	#plt.hist(a1, bins=bins)
	#plt.hist(a1, bins=np.arange(a1.min(), a1.max()-1))
	
	#plt.gca().set_xscale('log')


	#x2plot = np.arange(0, len(rmses), 1)
	#plt.plot(x2plot, rmses)
