#! /usr/bin/python

import scipy as sci
import scipy.special as sp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from math import acos,atan2

for N in range(1,100000,1000):
	Lmax = 20
	coeffs = np.zeros(Lmax + 1)

	vs = np.zeros((N,3))
	for s in range(0,N):
		while np.linalg.norm(vs[s], ord=2) > 1 or len(vs[s].nonzero()[0]) == 0:
			vs[s,0] = np.random.uniform(-1,1)
			vs[s,1] = np.random.uniform(-1,1)
			vs[s,2] = np.random.uniform(-1,1)

		vs[s] = vs[s]/np.linalg.norm(vs[s],ord=2)

	for l in range(0,Lmax + 1):
		ms = np.zeros(l*2+1)
		for m in range(-l,l+1):
			sumVal = 0
			count = 0
			for s in range(0,N):

				phi = np.arccos(vs[s,2])
				theta = np.arctan2(vs[s,1],vs[s,0]) + np.pi

				c = sp.sph_harm(m, l, theta, phi)
				sumVal += c
				count += 1

			ms[m+l] = sumVal*4*np.pi / count

		sumNorm = np.linalg.norm(ms,ord=2)
		coeffs[l] = sumNorm
		print 'l:' + str(l) + ' ' + str(coeffs[l])

	figure = plt.figure()
	plt.plot(np.arange(0,len(coeffs)),np.absolute(coeffs))
	plt.show()
