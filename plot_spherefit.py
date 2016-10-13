#! /usr/bin/python

import sys
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from src.OBJIO import loadOBJviaVTK
from opendr.topology import get_vert_connectivity
from pymetis import part_graph

if __name__ == '__main__':
	vertices, faces, normals, polyData = loadOBJviaVTK(sys.argv[1])

	cnct = get_vert_connectivity(vertices, faces)
	nbrs = [np.nonzero(np.array(cnct[:,i].todense()))[0] for i in range(cnct.shape[1])]

	edgeCuts, parts = part_graph(5,nbrs)

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')


	ax.scatter(vertices.T[0],vertices.T[1],vertices.T[2], color='b', marker='.')

	#ax.scatter(centerPoint[0],centerPoint[1],centerPoint[2], color='r',s=10)
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')
	plt.show()


