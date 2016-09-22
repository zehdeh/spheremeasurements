#! /usr/bin/python

import sys
import numpy as np
from src.shapes import Sphere
from opendr.topology import get_vert_connectivity
import matplotlib.pyplot as plt
from scipy.sparse import *
import scipy


if __name__ == '__main__':
	if len(sys.argv) < 2:
		print 'Please provide as an argument the file you want to plot'
		sys.exit(1)
	shape = Sphere(sys.argv[1],150)
	cnct = get_vert_connectivity(shape.vertices.T, shape.faces)
	nbrs = [np.nonzero(np.array(cnct[:,i].todense()))[0] for i in range(cnct.shape[1])]
	N = shape.vertices.shape[1]

	
	D = np.diag([len(d) for d in nbrs])

	A = cnct.toarray()
	A[A > 0] = 1
	L = D - A
	L_csc = bsr_matrix(L)

	k = 50#shape.vertices.shape[1]-20
	#eigenvals, eigenvectors = scipy.sparse.linalg.eigs(L_csc)
	eigenvals, eigenvectors = scipy.sparse.linalg.eigsh(L,k, which='SA')
	#eigenvals, eigenvectors = np.linalg.eig(L)
	eigenvals = np.real(eigenvals)
	eigenvectors = np.real(eigenvectors)

	#eigenvals,eigenvectors = np.linalg.eig(L)

	sortedEVindices = np.argsort(eigenvals)
	eigenvals = eigenvals[sortedEVindices]
	eigenvectors = eigenvectors[:,sortedEVindices]
	#print eigenvals

	Xhat = eigenvectors.T.dot(shape.vertices.T)

	indices = [i for i, f in enumerate(Xhat)]
	plt.plot(indices[:k], np.linalg.norm(Xhat,axis=1)[:k])
	plt.show()

	vertices2 = np.real(eigenvectors[:,:k].dot(Xhat[:k]))
	from mayavi import mlab

	mesh = mlab.triangular_mesh(vertices2.T[0], vertices2.T[1], vertices2.T[2], shape.faces)
	mlab.scalarbar(mesh)
	mlab.outline()
	mlab.show()
