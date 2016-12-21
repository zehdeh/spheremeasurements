USE_CACHE = False

import os
import numpy as np
import math
from scipy.special import sph_harm
from scipy.optimize import leastsq
from scipy.linalg import lstsq

def P(l, m, x):
	pmm = 1.0
	if m > 0:
		somx2 = np.sqrt((1.0 - x) * (1.0 + x))
		fact = 1.0
		for i in range(1, m+1):
			pmm *= (-fact)*somx2
			fact += 2.0

	if l == m:
		return pmm
	pmmp1 = x * (2.0 * m + 1.0)*pmm
	if l == (m + 1):
		return pmmp1

	pll = 0.0
	for ll in range(m+2, l+1):
		pll = ((2.0 * ll - 1.0) * x * pmmp1 - (ll + m - 1.0) * pmm) / (ll - m)
		pmm = pmmp1
		pmmp1 = pll
	
	return pll

def K(l, m):
	temp = ((2.0 * l + 1.0) * math.factorial(math.fabs(l - m))) / (4.0*np.pi*math.factorial(math.fabs(l+m)))
	return math.sqrt(math.fabs(temp))

def SH(l, m, theta, phi):
	sqrt2 = math.sqrt(2)
	if m == 0:
		return K(l, 0) * P(l, m, np.cos(theta))
	elif m > 0:
		return sqrt2*K(l,m) * np.cos(m*phi) * P(l, m, np.cos(theta))
	else:
		return sqrt2*K(l, -m) * np.sin(-m * phi) * P(l, -m, np.cos(theta))

def getSHMatrix(Lmax, phi, theta):
	Y_N = np.ones((phi.shape[0],(Lmax+1)**2),dtype=np.float)
	for l in range(Lmax + 1):
		Y_N[:,(l)**2:l**2+l*2+1] = np.array([SH(m-l, l, phi, theta).real for m in range(0, 2*l+1)]).T

	return Y_N

def back_transform(sphericalCoordinates, coefficients, Lmax, vertexAreas):
	phi, theta, radii = sphericalCoordinates

	reconstructedRadii = np.zeros((phi.shape[0]))
	totalArea = np.sum(vertexAreas)

	Y_N = getSHMatrix(Lmax, phi, theta)

	vals = list()
	#reconstructedRadii = Y_N.dot(coefficients)
	for i in range(phi.shape[0]):
		val = 0
		for l in range(Lmax):
			for m in range(-l, l+1):
				index = l*(l+1)+m
				val += coefficients[index] * Y_N[i,index]
		vals.append(val)
	reconstructedRadii = np.array(vals)

	return reconstructedRadii

def simple_transform(sphericalCoordinates, Lmax, vertexAreas):
	phi, theta, radii = sphericalCoordinates
	
	ys = np.zeros(Lmax + 1)
	coefficients = np.zeros((Lmax+1)**2)

	vertexAreas = (vertexAreas / np.sum(vertexAreas))*4*np.pi
	totalArea = np.sum(vertexAreas)

	#for l in range(Lmax):
	Y_N = getSHMatrix(Lmax, phi, theta)

	#coefficients = np.zeros(phi.shape[0])
	#for i in range((Lmax + 1)**2):
	#	print i
	#	coefficients[i] = np.sum([radii[j]*Y_N[j,i]*vertexAreas[j] for j in range(phi.shape[0])])

	for l in range(Lmax+1):
		for m in range(-l,l+1):
			coefficients[l*(l+1)+m] = np.sum([sph_harm(m,l, phi[i], theta[i]).real * r * vertexAreas[i] for i,r in enumerate(radii)])

	#coefficients = Y_N.T.dot(np.diag(vertexAreas)).dot(radii)
	print coefficients

	for l in range(Lmax + 1):
		ys[l] = np.linalg.norm(coefficients[l**2:(l**2 + 2*l + 1)], ord=2)
	return ys, coefficients

def matlab_approach(sphericalCoordinates, Lmax, vertexAreas):
	phi, theta, radii = sphericalCoordinates
	#radii = np.ones(theta.shape[0])
	Y_N = np.ones((phi.shape[0],(Lmax+1)**2),dtype=np.float)
	for l in range(Lmax + 1):
		Y_N[:,(l)**2:l**2+l*2+1] = np.array([sph_harm(m-l, l, phi, theta) for m in range(0, 2*l+1)]).T
	
	#b = np.linalg.inv(Y_N.T.dot(Y_N)).dot(Y_N.T).dot(radii)
	b = lstsq(Y_N*np.sqrt(vertexAreas)[:,None],radii*np.sqrt(vertexAreas))[0]
	return Y_N,b

def residuals(aa, Yval, rs):
	return (Yval.dot(aa) - rs)

def getSphericalCoordinates(vertices, centerPoint):
	vertices = vertices - centerPoint

	#vertices = vertices / np.linalg.norm(vertices, axis=1,ord=2)[...,None]

	r = np.sqrt(vertices.T[0]**2 + vertices.T[1]**2 + vertices.T[2]**2)
	theta = np.arccos(vertices.T[2]/r)
	phi = np.arctan2(vertices.T[1],vertices.T[0]) + np.pi
	print r.min(), r.max()
	return np.array([phi, theta, r])

def getCartesianCoordinates(theta, phi, r, centerPoint):
	x = r*np.sin(theta)*np.cos(phi)
	y = r*np.sin(theta)*np.sin(phi)
	z = r*np.cos(theta)

	coords = np.array([x,y,z])
	coords = coords.T + centerPoint

	return coords

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
				s = sph_harm(l, m, phi[i], theta[i]).real
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

def getAMatrix(theta, phi, l, filePath, oldA=None):
	n = len(phi)
	
	path = os.path.split(filePath)
	fileName = path[-1][:-4]

	matrixName = fileName + '_' + str(l)
	matrixPath = os.path.join(path[0],matrixName)
	if os.path.isfile(matrixPath + '.npy') and USE_CACHE:
		A = np.load(matrixPath + '.npy')
		print 'Using cache for ' + str(l)
	else:
		A = np.zeros((n, 2*l+1), dtype=np.complex)
		if oldA is None:
			for i in range(n):
				A[i] = [sph_harm(j-l, l, phi[i], theta[i]) for j in range(0, 2*l+1)]
		else:
			A[:,1:2*l] = oldA
			for i in range(n):
				A[i,0] = sph_harm(-l, l, phi[i], theta[i])
				A[i,-1] = sph_harm(l, l, phi[i], theta[i])
		np.save(matrixPath,A)
	
	return A
def sirf(sphericalCoordinates, Lmax, fileName, vertexAreas):
	phi, theta, r = sphericalCoordinates

	s = min(5, Lmax)
	# areas auf beide Seiten!
	A = getAMatrix(phi, theta, 0, fileName)
	for i in range(1,s):
		A = np.concatenate((A,getAMatrix(phi, theta, i, fileName)), axis=1)
	
	A = A*np.sqrt(vertexAreas)[:,None]
	b = np.zeros((Lmax+1)**2, dtype=np.complex)
	#b[0:A.shape[1]], flag = leastsq(residuals, b[0:A.shape[1]], args=(A,r))
	b[0:A.shape[1]] = lstsq(A,r)[0]

	np.set_printoptions(threshold=np.inf)
	res = r - A.dot(b[0:A.shape[1]])
	Al = getAMatrix(phi, theta, s-1, fileName)
	for l in range(s, Lmax+1):
		Al = getAMatrix(phi, theta, l, fileName, Al)
		lowerBound = l**2
		upperBound = (l+1)**2

		#b[lowerBound:upperBound], flag = leastsq(residuals, b[lowerBound:upperBound], args=(Al,res))
		b[lowerBound:upperBound] = lstsq(Al,res)[0]
		res = res - Al.dot(b[lowerBound:upperBound])
		A = np.concatenate((A,Al), axis=1)

	ys = np.zeros(Lmax+1)
	for i in range(Lmax+1):
		lowerBound = i**2
		upperBound = (i+1)**2
		ys[i] = np.linalg.norm(b[lowerBound:upperBound], ord=2)
	
	return b, ys, A
