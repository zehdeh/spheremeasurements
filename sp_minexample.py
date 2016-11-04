#! /usr/bin/python

import scipy as sci
import scipy.special as sp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from math import acos,atan2
from src.OBJIO import loadOBJ, getBounds
from src.fitting import fitSphere

Lmax = 20
coeffs = np.zeros(Lmax + 1)

vertices, faces, normals = loadOBJ('res/fine/synthetic_150_fine.obj')
bounds = getBounds(vertices)
centerPoint, radius = fitSphere(vertices, [0,0,0,1],150, bounds)
vertices = vertices - centerPoint
#vertices = vertices / vertices.max()
vertices = vertices / np.linalg.norm(vertices, axis=1,ord=2)[...,None]

N = 3000
vs = np.zeros((N,3))
for s in range(0,N):
	while np.linalg.norm(vs[s], ord=2) > 1 or len(vs[s].nonzero()[0]) == 0:
		vs[s,0] = np.random.uniform(-1,1)
		vs[s,1] = np.random.uniform(-1,1)
		vs[s,2] = np.random.uniform(-1,1)

	vs[s] = vs[s]/np.linalg.norm(vs[s],ord=2)

vs = vertices

phi = np.arccos(vs[:,2])
theta = np.arctan2(vs[:,1],vs[:,0]) + np.pi
for l in range(0,Lmax + 1):
	ms = [np.sum(sp.sph_harm(m, l, theta, phi))*4*np.pi / vs.shape[0] for m in range(-l,l+1)]

	sumNorm = np.linalg.norm(ms,ord=2)
	coeffs[l] = sumNorm
	print 'l:' + str(l) + ' ' + str(coeffs[l])

figure = plt.figure()
plt.plot(np.arange(0,len(coeffs)),np.absolute(coeffs))
plt.show()
