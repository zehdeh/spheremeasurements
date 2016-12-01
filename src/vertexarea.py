import numpy as np
import math

def getFaceArea(face, vertices):
	v = vertices[face]
	a = np.linalg.norm(v[0] - v[1])
	b = np.linalg.norm(v[1] - v[2])
	c = np.linalg.norm(v[2] - v[0])

	t = [a,b,c]

	s = np.sum(t)/2
	return np.sqrt(s*np.prod(s-t))

def getVertexAreas(faces, vertices):
	faceAreas = np.array([getFaceArea(f, vertices) for f in faces])

	vertexAreas = np.array([np.sum(faceAreas[np.argwhere(np.any(np.reshape(i == faces.ravel(), (-1,3)), axis=1))]/3) for i,v in enumerate(vertices)])
	
	return vertexAreas
