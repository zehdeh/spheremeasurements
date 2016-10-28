import numpy as np
from math import acos
def getVoronoiArea(vertices, faces):
	voronoiVertices = np.zeros((faces.shape[0],3))
	for n in range(faces.shape[0]):
		i0 = faces[n, 0]
		i1 = faces[n, 1]
		i2 = faces[n, 2]

		r_01 = vertices[i1] - vertices[i0]
		r_02 = vertices[i2] - vertices[i0]

		r_normal = np.cross(r_01, r_02)
		u_normal = r_normal/np.linalg.norm(r_normal)

		voronoiVertices[n] = u_normal
	
	for n in range(vertices.shape[0]):
		faceIdx = np.where(faces == n)

		k = 0
		currentFaceIdx = faceIdx[k]
		currentFace = faces[currentFaceIdx]
		currentVertexIdx = np.where(currentFace != n)
		currentVertex = currentFace[currentVertexIdx]

		sortedFaceIdx = faceIdx[k]
		notSorted = True

		while notSorted:
			tempFaceList = faceIdx(np.where(currentFaceIdx != faceIdx))

			for l in range(len(tempFaceList)):
				currentFaceIdx = tempFaceList[l]
				currentFace = faces[currentFaceIdx]

				if np.any(currentFace == currentVertex):
					sortedFaceIdx.append(currentFaceIdx)
					if len(sortedFaceIdx) == len(faceIdx):
						notSorted = False
						break

					currentVertexIdx = np.where((currentFace != n)*(currentFace != currentVertex))
					currentVertex = currentFace[currentVertexIdx]

		for i in range(len(sortedFaceIdx)):
			if duplicates[sortedFaceIdx[i]] != 0:
				sortedFaceIdx[i] = duplicates[sortedFaceIdx[i]]

	
	areas = np.zeros(voronoiFaces.shape[0])
	for m,face in enumerate(voronoiFaces):
		theta = np.zeros(voronoiVertices.shape[0])
		for n in range(3):
			r_01 = voronoiVertices[face[0]]
			r_02 = voronoiVertices[face[1]]

			r_2x1 = np.cross(r_02, r_01)

			r_21 = np.cross(r_2x1, r_02)

			r_03 = voronoiVertices[face[2]]

			r_2x3 = np.cross(r_02, r_03)

			r_23 = np.cross(r_2x3, r_02)

			r_21 = r_21/np.linalg.norm(r_21)
			r_23 = r_23/np.linalg.norm(r_23)
			theta[n] = acos(np.dot(r_21,r_23))

			face = np.roll(face, -1)

		areas[m] = np.sum(theta) - np.pi
	return areas
