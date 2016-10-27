import numpy as np
from math import acos
def getVoronoiArea(vertices, faces):
	areas = np.zeros(faces.shape[0])
	for m,face in enumerate(faces):
		theta = np.zeros(vertices.shape[0])
		for n in range(3):
			r_01 = vertices[face[0]]
			r_02 = vertices[face[1]]

			r_2x1 = np.cross(r_02, r_01)

			r_21 = np.cross(r_2x1, r_02)

			r_03 = vertices[face[2]]

			r_2x3 = np.cross(r_02, r_03)

			r_23 = np.cross(r_2x3, r_02)

			r_21 = r_21/np.linalg.norm(r_21)
			r_23 = r_23/np.linalg.norm(r_23)
			theta[n] = acos(np.dot(r_21,r_23))

			face = np.roll(face, -1)

		areas[m] = np.sum(theta) - np.pi
	return areas
